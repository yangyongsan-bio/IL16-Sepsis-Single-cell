#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
8.0.prepare_validation_inputs.py
================================================================================
独立验证队列（GSE65682 成人 / GSE66099 儿童）
—— GEO 原始数据下载 + 8.1 / 8.2 输入文件构建 ——

【用途】
把「从 GEO 下载」到「8.1 / 8.2 可直接运行」之间的全部中间步骤脚本化，
使独立验证队列的流程完全可复现。

【产出：即 8.1 / 8.2 的输入】
  GSE65682_validation/
    GSE65682_series_matrix.txt.gz      ← 下载（GEO 原始；其表体已由作者完成
                                          RMA + 0.5 方差过滤(49386→24646) + ComBat，
                                          见文件头 !Sample_data_processing）
    GSE65682_expr_matrix.csv.gz        ← 构建（24646 探针 × 802 样本）
    GSE65682_pheno.tsv                 ← 构建（802 样本 × 11 列）
    gene_probe_map.json                ← 构建（17 基因 → GPL13667/U219 探针）
  GSE66099_validation/
    GSE66099_series_matrix.txt.gz      ← 下载（8.2 直接读它，不做转换）
    subseries/GSE26378_series_matrix.txt.gz   ← 下载（结局标签来源）
    subseries/GSE26440_series_matrix.txt.gz   ← 下载（结局标签来源）
    gsm_to_original.json               ← 构建（276 样本 "Reanalysis of" 映射）
    pediatric_outcome_map.json         ← 构建（子系列 28 天结局标签）
    gene_probe_map_ped.json            ← 构建（21 基因 → GPL570/U133Plus2 探针）
================================================================================
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import socket
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime

# Windows 控制台中文输出
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = os.path.dirname(os.path.abspath(__file__))
AD = os.path.join(ROOT, "GSE65682_validation")     # 成人
PD = os.path.join(ROOT, "GSE66099_validation")     # 儿童
SUBSERIES = os.path.join(PD, "subseries")

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"

# ---- 下载清单: GSE -> 本地目标路径 ----
DOWNLOADS = {
    "GSE65682": os.path.join(AD, "GSE65682_series_matrix.txt.gz"),
    "GSE66099": os.path.join(PD, "GSE66099_series_matrix.txt.gz"),
    "GSE26378": os.path.join(SUBSERIES, "GSE26378_series_matrix.txt.gz"),
    "GSE26440": os.path.join(SUBSERIES, "GSE26440_series_matrix.txt.gz"),
}
# 注：GSE66099 是 SuperSeries（其 !Series_summary 点名汇总了 6 个子库），
#     但只有 GSE26378 / GSE26440 带显式 outcome: 字段（= 正文 108 例的来源）

# ---- 注释包 ----
HG_U219 = os.path.join(AD, "hgu219.db", "inst", "extdata", "hgu219.sqlite")
HG_U133P2 = os.path.join(PD, "hgu133plus2.db", "inst", "extdata", "hgu133plus2.sqlite")
ADULT_ENTREZ = {
    # 成人 T CD8 TEMRA
    "UPF3A": 65110, "SOCS1": 8651, "CYB561D1": 284613, "PLIN2": 123,
    # 成人髓系 cDC2
    "MGA": 23269, "MTMR9": 66036, "SOD1": 6647,
    # 通讯配体
    "GAS6": 2621, "IL16": 3603,
    # 跨模块参考基因（SMR / 儿童模块）
    "CBX1": 10951, "FBXO7": 25793, "TNIP2": 79155,
    "UFD1": 7353, "IMMT": 10989, "THOP1": 7064,
    # 预期无探针
    "ND4L": 4539, "LINC00211": None,
}

PED_ENTREZ = {
    # 儿童髓系 Inflammatory Monocytes（16）
    "WNK1": 65125, "CCT7": 10574, "HDGF": 3068, "FBXO7": 25793, "TNIP2": 79155,
    "UBE3C": 9690, "CBX1": 10951, "TRPV2": 51393, "VGLL4": 9686, "SUN2": 25777,
    "MVD": 4597, "IARS2": 55699, "U2AF2": 11338, "DCXR": 51181, "CSNK2A2": 1459,
    "MPHOSPH6": 10200,
    # 儿童 T CD4 Early Activated（3）
    "UFD1": 7353, "IMMT": 10989, "THOP1": 7064,
    # 跨谱系 hub + lncRNA
    "IL16": 3603, "LINC00211": None,
}

EXPECT = dict(
    adult_probes=24646, adult_samples=802, adult_outcome_n=479, adult_death=114,
    adult_assessable=9,                    # 11 目标基因中，9 个有探针（ND4L / LINC00211 无）
    ped_samples=276, ped_shock=181, ped_matched=108, ped_death=13, ped_assessable=20,
    sub_gse26378=103, sub_gse26440=130,
)


ADULT_TARGET = ["UPF3A", "SOCS1", "CYB561D1", "PLIN2", "ND4L", "MGA", "MTMR9",
                "SOD1", "IL16", "GAS6", "LINC00211"]
PED_TARGET = list(PED_ENTREZ.keys())        # 8.2 遍历 JSON 的全部 21 个键
NO_PROBE = {"ND4L", "LINC00211"}            # 平台无探针 → not assessable（预期）

BACKUP_DIR = None    # 由 setup_backup() 填充
FORCE_WRITE = False  # 由 --force-write 置 True：无条件重写产出

# 1. 通用工具
def log(msg=""):
    print(msg, flush=True)


def hr(title=""):
    log("\n" + "=" * 80)
    if title:
        log("■ " + title)
        log("=" * 80)


def geo_prefix(gse):
    """GSE 号 -> GEO FTP 目录前缀，如 GSE65682 -> GSE65nnn, GSE4607 -> GSE4nnn,
    GSE1000 -> GSEnnn"""
    d = gse[3:]
    return "GSE" + d[:-3] + "nnn" if len(d) > 3 else "GSEnnn"


def http_size(url, tries=3):
    """HEAD 请求取远端大小；失败返回 None。"""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, method="HEAD",
                                         headers={"User-Agent": "Mozilla/5.0"})
            r = urllib.request.urlopen(req, timeout=30)
            return int(r.headers.get("Content-Length") or 0), r.geturl()
        except Exception as e:
            if k == tries - 1:
                log("    HEAD 失败: %s: %s" % (type(e).__name__, e))
                return None, None
            time.sleep(2 * (k + 1))


def download(gse, dest, force=False):
    """下载一个 GSE 的 series matrix。返回 True/False。"""
    if os.path.exists(dest) and not force:
        size = os.path.getsize(dest)
        url = GEO.format(prefix=geo_prefix(gse), gse=gse)
        rsize, real = http_size(url)
        if rsize and rsize == size:
            log("  %-9s 已存在且大小一致 (%d B)，跳过" % (gse, size))
            return True
        if rsize:
            log("  %-9s 已存在但大小不符（本地 %d / 远端 %d）→ 重新下载" % (gse, size, rsize))
        else:
            log("  %-9s 已存在（远端大小未知），跳过" % gse)
            return True

    cands = [GEO.format(prefix=geo_prefix(gse), gse=gse)]
    alt = GEO.format(prefix="GSEnnn", gse=gse)
    if alt not in cands:
        cands.append(alt)

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for url in cands:
        try:
            log("  下载 %s ..." % url)
            t0 = time.time()
            with urllib.request.urlopen(
                    urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}),
                    timeout=180) as r, open(dest + ".part", "wb") as f:
                shutil.copyfileobj(r, f, 1024 * 1024)
            os.replace(dest + ".part", dest)
            log("    OK  %.1f MB  用时 %.0fs" % (os.path.getsize(dest) / 1048576, time.time() - t0))
            return True
        except Exception as e:
            log("    ✗ %s: %s" % (type(e).__name__, e))
            if os.path.exists(dest + ".part"):
                os.remove(dest + ".part")
    return False


# ---- series matrix 解析 ----
def gz_lines(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace")


def parse_header(path):
    """返回 (field_lines, all_lines)
    field_lines: list[(field_name, [value per sample])]   —— 按出现顺序保留重复字段
    all_lines  : 头部所有行
    """
    fields, raw = [], []
    with gz_lines(path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith("!series_matrix_table_begin"):
                break
            raw.append(line)
            if line.startswith("!Sample_"):
                parts = line.split("\t")
                fields.append((parts[0][len("!Sample_"):],
                               [v.strip('"') for v in parts[1:]]))
    return fields, raw


def field_map(path):
    """{字段名: [每样本值]}；同名多行自动加 _2 / _3 后缀（characteristics 尤其注意）。"""
    fields, _ = parse_header(path)
    out, seen = {}, {}
    for name, vals in fields:
        seen[name] = seen.get(name, 0) + 1
        key = name if seen[name] == 1 else "%s_%d" % (name, seen[name])
        out[key] = vals
    return out


def char_dict(path):
    """把所有 !Sample_characteristics_ch1 行拆成 {字段: {GSM: 值}}"""
    fields, _ = parse_header(path)
    gsm = None
    for name, vals in fields:
        if name == "geo_accession":
            gsm = vals
            break
    if gsm is None:
        raise RuntimeError("找不到 !Sample_geo_accession")
    res = {}
    for name, vals in fields:
        if not name.startswith("characteristics"):
            continue
        for g, raw in zip(gsm, vals):
            raw = raw.strip()
            if ":" not in raw:
                continue
            k, v = raw.split(":", 1)
            res.setdefault(k.strip(), {})[g] = v.strip()
    return gsm, res


def table_meta(path):
    """表体维度 + 探针 ID 列表（流式，省内存）。返回 (probe_ids, sample_names)"""
    probes, samples = [], []
    with gz_lines(path) as f:
        started = False
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if not started:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not line.strip():
                continue
            parts = line.split("\t")
            if parts[0].strip('"') == "ID_REF":
                samples = [p.strip('"') for p in parts[1:]]
                continue
            probes.append(parts[0].strip('"'))
    return probes, samples


# ---- 内容级比对与写盘 ----
def digest(path, gz=False):
    h = hashlib.md5()
    if gz:
        with gzip.open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    else:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    return h.hexdigest()


def same_content(a, b, gz=False):
    """比较两个文件的内容是否一致。
    · gz=True：比解压后的字节流（gzip 容器含 mtime，逐字节不可比）
    · 先比字节；不同则忽略行尾差异（LF vs CRLF）再比一次
      （pandas / json 在不同版本与平台下行尾不同 —— 本机历史产物是 CRLF，
        Python 文本/二进制写入是 LF；对 R 的 fread 无任何影响）
    返回 (内容一致, 字节完全一致)
    """
    def reader(p):
        return gzip.open(p, "rb") if gz else open(p, "rb")
    with reader(a) as fa, reader(b) as fb:
        ba, bb = fa.read(), fb.read()
    if ba == bb:
        return True, True
    return ba.replace(b"\r\n", b"\n") == bb.replace(b"\r\n", b"\n"), False


def commit(tmp, dest, gz=False, label="", force=None):
    """落盘。返回 (status, detail)。
    · force=True：无条件重写（仍先备份）—— 纯「生成」模式
    · force=False：与现有一致则不动（保留原文件时间戳），不一致才备份覆盖
    · force=None（默认）：取全局 FORCE_WRITE（--force-write 开关）
    """
    global BACKUP_DIR
    if force is None:
        force = FORCE_WRITE
    if os.path.exists(dest):
        same, ident = (False, False) if force else same_content(tmp, dest, gz=gz)
        if same:
            os.remove(tmp)
            return ("same", "与现有文件内容一致，未改动" if ident
                    else "内容一致（行尾 CRLF/LF 差异已忽略），未改动")
        if BACKUP_DIR:
            os.makedirs(BACKUP_DIR, exist_ok=True)
            shutil.copy2(dest, os.path.join(BACKUP_DIR, os.path.basename(dest)))
        os.replace(tmp, dest)
        return "updated", "内容有变化 → 已备份并覆盖"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    os.replace(tmp, dest)
    return "created", "新建"


def report(label, status, detail, extra=""):
    icon = {"same": "✅", "updated": "🔁", "created": "🆕"}[status]
    log("  %s %-34s %s%s" % (icon, label, detail, ("  " + extra) if extra else ""))


# 2. 成人（GSE65682）：expr_matrix.csv.gz + pheno.tsv
def build_adult_expr(src_gz, dest, tmp):
    n_probe = 0
    with gz_lines(src_gz) as fin, gzip.open(tmp, "wb") as fout:
        started = False
        for line in fin:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if not started:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not line.strip():
                continue
            parts = [p.strip('"') for p in line.split("\t")]
            if parts[0] == "ID_REF":
                parts[0] = ""            # 与 pandas 产物一致：索引列无列名
            else:
                n_probe += 1
            fout.write((",".join(parts) + "\n").encode("utf-8"))
    return n_probe


def build_adult_pheno(src_gz, dest, tmp):
    
    import numpy as np
    import pandas as pd

    gsm_order, chars = char_dict(src_gz)
    pheno = pd.DataFrame({"gsm": gsm_order}).set_index("gsm")

    def col(key):
        return pd.Series(pheno.index.map(chars.get(key, {})), index=pheno.index)

    pheno["gender"] = col("gender")
    pheno["age"] = pd.to_numeric(col("age").replace({"NA": np.nan}), errors="coerce")
    pheno["pneumonia_dx"] = col("pneumonia diagnoses")
    pheno["mortality_28"] = pd.to_numeric(
        col("mortality_event_28days").replace({"NA": np.nan}), errors="coerce")
    pheno["ttd_28"] = pd.to_numeric(
        col("time_to_event_28days").replace({"NA": np.nan}), errors="coerce")
    pheno["endotype_cohort"] = col("endotype_cohort")
    pheno["endotype_class"] = col("endotype_class")
    pheno["icu_infection"] = col("icu_acquired_infection")
    pheno["diabetes"] = col("diabetes_mellitus")
    pheno["abdominal"] = col("abdominal_sepsis_and_controls")

    pheno.to_csv(tmp, sep="\t")
    return pheno

# 3. 探针注释
def entrez_to_probes(sqlite_path):
    if not os.path.exists(sqlite_path):
        raise FileNotFoundError(
            "缺少注释包: %s\n"
            "请确认验证目录下存在 Bioconductor 注释包（hgu219.db / hgu133plus2.db）"
            % sqlite_path)
    con = sqlite3.connect(sqlite_path)
    by_gid = {}
    for pid, gid in con.execute("SELECT probe_id, gene_id FROM probes"):
        by_gid.setdefault(str(gid), []).append(pid)
    con.close()
    return by_gid


def build_probe_map(source, sqlite_path, dest, tmp, entrez_as_str=False):
    """构建 {基因: {entrez, probes}}。entrez_as_str=True 时 Entrez 写为字符串（与儿童历史产物一致）。"""
    by_gid = entrez_to_probes(sqlite_path)
    out = {}
    for sym, eg in source.items():
        probes = by_gid.get(str(eg), []) if eg is not None else []
        out[sym] = {"entrez": (str(eg) if (eg is not None and entrez_as_str) else eg),
                    "probes": probes}
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out


def refresh_entrez(symbols):
    """用 mygene.info 复核 symbol -> Entrez（可选）。返回 {symbol: entrez or None}"""
    import urllib.parse
    q = urllib.parse.quote(" OR ".join("symbol:%s" % s for s in symbols))
    url = ("https://mygene.info/v3/query?q=%s&species=human"
           "&fields=symbol,entrezgene&size=100" % q)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    d = json.loads(urllib.request.urlopen(req, timeout=60).read())
    got = {}
    for h in d.get("hits", []):
        if h.get("symbol") and h.get("entrezgene"):
            got[h["symbol"]] = int(h["entrezgene"])
    return got


# 4. 儿童（GSE66099）：gsm_to_original.json + pediatric_outcome_map.json
def build_gsm_map(src_gz, dest, tmp):
    fields, _ = parse_header(src_gz)
    gsm = None
    for name, vals in fields:
        if name == "geo_accession":
            gsm = vals
            break
    disease, original = {}, {}
    for name, vals in fields:
        if name.startswith("characteristics"):
            for g, raw in zip(gsm, vals):
                if raw.startswith("disease:"):
                    disease[g] = raw.split(":", 1)[1].strip()
        elif name.startswith("relation"):
            for g, raw in zip(gsm, vals):
                m = re.search(r"Reanalysis of:\s*(GSM\d+)", raw)
                if m:
                    original[g] = m.group(1)
    out = {"gsm_order": gsm, "disease": disease, "original": original}
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(out, f)          # 单行，与历史产物一致
    return out


def build_ped_outcome(sub_files, dest, tmp):
    out = {}
    for gse in sub_files:
        p = DOWNLOADS.get(gse) or os.path.join(SUBSERIES, "%s_series_matrix.txt.gz" % gse)
        if not os.path.exists(p):
            raise FileNotFoundError("缺少子系列原始文件: %s" % p)
        gsm, chars = char_dict(p)
        oc = chars.get("outcome", {})
        out[gse] = {g: oc[g] for g in gsm if g in oc}
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)   # 与历史产物一致
    return out


# 5. 自检
def verify():
    hr("一致性自检")
    ok = True

    def chk(cond, label, got, want):
        nonlocal ok
        mark = "✅" if cond else "❌"
        if not cond:
            ok = False
        log("  %s %-46s 实测 %-10s 期望 %s" % (mark, label, got, want))

    # ---- 成人 ----
    log("\n[成人 GSE65682]")
    ad_gz = DOWNLOADS["GSE65682"]
    if os.path.exists(ad_gz):
        probes, samples = table_meta(ad_gz)
        chk(len(probes) == EXPECT["adult_probes"], "series matrix 探针数", len(probes), EXPECT["adult_probes"])
        chk(len(samples) == EXPECT["adult_samples"], "series matrix 样本数", len(samples), EXPECT["adult_samples"])
    ex = os.path.join(AD, "GSE65682_expr_matrix.csv.gz")
    ph = os.path.join(AD, "GSE65682_pheno.tsv")
    if os.path.exists(ex) and os.path.exists(ph):
        import csv as _csv
        with gzip.open(ex, "rt", encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split(",")
            n_expr = sum(1 for _ in f)
        chk(n_expr == EXPECT["adult_probes"], "expr_matrix 探针数", n_expr, EXPECT["adult_probes"])
        chk(len(hdr) - 1 == EXPECT["adult_samples"], "expr_matrix 样本数", len(hdr) - 1, EXPECT["adult_samples"])
        rows = list(_csv.reader(open(ph, encoding="utf-8"), delimiter="\t"))
        cols = rows[0]
        idx = {c: i for i, c in enumerate(cols)}
        m28 = [r[idx["mortality_28"]] for r in rows[1:]]
        n_out = sum(1 for v in m28 if v != "")
        n_dead = sum(1 for v in m28 if v.startswith("1"))
        chk(len(rows) - 1 == EXPECT["adult_samples"], "pheno 样本数", len(rows) - 1, EXPECT["adult_samples"])
        chk(n_out == EXPECT["adult_outcome_n"], "pheno 有 28 天结局的样本", n_out, EXPECT["adult_outcome_n"])
        chk(n_dead == EXPECT["adult_death"], "pheno 28 天死亡数", n_dead, EXPECT["adult_death"])
        chk([r[idx["gsm"]] for r in rows[1:]] == hdr[1:], "pheno 样本顺序 == expr 列顺序", "一致", "一致")
    gm = os.path.join(AD, "gene_probe_map.json")
    if os.path.exists(gm) and os.path.exists(ex):
        d = json.load(open(gm, encoding="utf-8"))
        have = set(hdr[1:]) if False else set()
        with gzip.open(ex, "rt", encoding="utf-8") as f:
            f.readline()
            have = {ln.split(",", 1)[0] for ln in f}
        log("\n  探针覆盖（平台注释 vs 方差过滤后的矩阵；矩阵 = 24646 探针）：")
        log("  %-11s %-8s %-8s %s" % ("gene", "平台", "矩阵内", "状态"))
        n_assess = 0
        bad = {}
        for g in ADULT_TARGET:
            ps = d.get(g, {}).get("probes", [])
            inn = [p for p in ps if p in have]
            if not ps:
                tag = "不可评估（平台无探针）"
            elif inn:
                tag = "可评估"; n_assess += 1
            else:
                tag = "❌ 平台有探针但全被方差过滤"
                bad[g] = ps
            log("  %-11s %-8d %-8d %s" % (g, len(ps), len(inn), tag))
        chk(not bad, "8.1 的 11 个目标基因均可评估或确证无探针",
            "全部通过" if not bad else bad, "全部通过")
        chk(n_assess == EXPECT["adult_assessable"], "成人可评估基因数",
            n_assess, EXPECT["adult_assessable"])

    # ---- 儿童 ----
    log("\n[儿童 GSE66099]")
    pd_gz = DOWNLOADS["GSE66099"]
    if os.path.exists(pd_gz):
        probes9, samples9 = table_meta(pd_gz)
        chk(len(samples9) == EXPECT["ped_samples"], "series matrix 样本数", len(samples9), EXPECT["ped_samples"])
        chk(len(probes9) == 54675, "series matrix 探针数", len(probes9), 54675)
    g2o = os.path.join(PD, "gsm_to_original.json")
    om = os.path.join(PD, "pediatric_outcome_map.json")
    if os.path.exists(g2o) and os.path.exists(om):
        g = json.load(open(g2o, encoding="utf-8"))
        o = json.load(open(om, encoding="utf-8"))
        chk(len(g["gsm_order"]) == EXPECT["ped_samples"], "gsm_to_original 样本数",
            len(g["gsm_order"]), EXPECT["ped_samples"])
        shock = [x for x in g["gsm_order"]
                 if g["disease"].get(x) == "SepticShock" and x in g["original"]]
        chk(len(shock) == EXPECT["ped_shock"], "SepticShock 样本数", len(shock), EXPECT["ped_shock"])
        matched = {}
        for x in shock:
            org = g["original"][x]
            for gse in o:                       # 与 8.2 相同遍历顺序（GSE26378 先）
                if org in o[gse]:
                    v = o[gse][org]
                    if v in ("Survivor", "Nonsurvivor"):
                        matched[x] = 1 if v == "Nonsurvivor" else 0
                    break
        chk(len(matched) == EXPECT["ped_matched"], "匹配到结局的样本", len(matched), EXPECT["ped_matched"])
        chk(sum(matched.values()) == EXPECT["ped_death"], "其中死亡数",
            sum(matched.values()), EXPECT["ped_death"])
        chk(len(o.get("GSE26378", {})) == EXPECT["sub_gse26378"], "GSE26378 结局标签数",
            len(o.get("GSE26378", {})), EXPECT["sub_gse26378"])
        chk(len(o.get("GSE26440", {})) == EXPECT["sub_gse26440"], "GSE26440 结局标签数",
            len(o.get("GSE26440", {})), EXPECT["sub_gse26440"])
    gmp = os.path.join(PD, "gene_probe_map_ped.json")
    if os.path.exists(gmp) and os.path.exists(pd_gz):
        d9 = json.load(open(gmp, encoding="utf-8"))
        have9 = set(probes9)
        log("\n  探针覆盖（GPL570 注释 vs gcRMA 矩阵；矩阵 = 54675 探针，未过滤）：")
        log("  %-11s %-8s %-8s %s" % ("gene", "平台", "矩阵内", "状态"))
        n_assess9 = 0
        bad9 = {}
        for g in PED_TARGET:
            ps = d9.get(g, {}).get("probes", [])
            inn = [p for p in ps if p in have9]
            if not ps:
                tag = "不可评估（平台无探针）"
            elif inn:
                tag = "可评估"; n_assess9 += 1
            else:
                tag = "❌ 平台有探针但矩阵内没有"
                bad9[g] = ps
            log("  %-11s %-8d %-8d %s" % (g, len(ps), len(inn), tag))
        chk(not bad9, "8.2 的 21 个目标基因均可评估或确证无探针",
            "全部通过" if not bad9 else bad9, "全部通过")
        chk(n_assess9 == EXPECT["ped_assessable"], "儿童可评估基因数",
            n_assess9, EXPECT["ped_assessable"])

    log("")
    log("  总结: %s" % ("✅ 全部通过" if ok else "❌ 存在不通过项，请检查上方 ❌"))
    return ok

# 6. 主流程
def setup_backup(enable):
    global BACKUP_DIR
    if enable:
        BACKUP_DIR = os.path.join(ROOT, "_待清理",
                                 "9.0重建_%s" % datetime.now().strftime("%Y%m%d_%H%M"))
    else:
        BACKUP_DIR = None


def do_download(args):
    hr("① 下载 GEO 原始数据")
    items = dict(DOWNLOADS)
    bad = []
    for gse, dest in items.items():
        if not download(gse, dest, force=args.force):
            bad.append(gse)
    if bad:
        log("\n  ❌ 下载失败: %s" % ", ".join(bad))
        log("     若为本机网络/代理问题，可重跑；或手动从 %s 下载后放到对应路径。"
            % "https://www.ncbi.nlm.nih.gov/geo/")
        return False
    return True


def do_build(args, tmpdir):
    hr("② 构建 8.1 / 8.2 的输入文件")
    results = []

    # ---------- 成人 ----------
    log("\n[成人 GSE65682]")
    src = DOWNLOADS["GSE65682"]
    if not os.path.exists(src):
        log("  ❌ 缺少原始文件 %s，请先运行 download" % src)
        return False

    tmp = os.path.join(tmpdir, "expr.csv.gz")
    n = build_adult_expr(src, None, tmp)
    st, dt = commit(tmp, os.path.join(AD, "GSE65682_expr_matrix.csv.gz"), gz=True)
    report("GSE65682_expr_matrix.csv.gz", st, dt, "(%d 探针)" % n)
    results.append(("expr", n))

    tmp = os.path.join(tmpdir, "pheno.tsv")
    pheno = build_adult_pheno(src, None, tmp)
    st, dt = commit(tmp, os.path.join(AD, "GSE65682_pheno.tsv"))
    report("GSE65682_pheno.tsv", st, dt, "(%d 样本)" % pheno.shape[0])
    results.append(("pheno", pheno.shape[0]))

    tmp = os.path.join(tmpdir, "gpm_adult.json")
    d = build_probe_map(ADULT_ENTREZ, HG_U219, None, tmp, entrez_as_str=False)
    st, dt = commit(tmp, os.path.join(AD, "gene_probe_map.json"))
    report("gene_probe_map.json", st, dt,
           "(%d 基因, 有探针 %d)" % (len(d), sum(1 for v in d.values() if v["probes"])))

    # ---------- 儿童 ----------
    log("\n[儿童 GSE66099]")
    src9 = DOWNLOADS["GSE66099"]
    if not os.path.exists(src9):
        log("  ❌ 缺少原始文件 %s，请先运行 download" % src9)
        return False
    log("  ✅ %-34s 已就位（8.2 直接读该文件，无需转换）" % "GSE66099_series_matrix.txt.gz")

    tmp = os.path.join(tmpdir, "gsm_to_original.json")
    g = build_gsm_map(src9, None, tmp)
    shock = [x for x in g["gsm_order"] if g["disease"].get(x) == "SepticShock"]
    st, dt = commit(tmp, os.path.join(PD, "gsm_to_original.json"))
    report("gsm_to_original.json", st, dt,
           "(%d 样本 / %d SepticShock)" % (len(g["gsm_order"]), len(shock)))

    tmp = os.path.join(tmpdir, "outcome.json")
    o = build_ped_outcome(["GSE26378", "GSE26440"], None, tmp)
    st, dt = commit(tmp, os.path.join(PD, "pediatric_outcome_map.json"))
    report("pediatric_outcome_map.json", st, dt,
           "(%s)" % ", ".join("%s %d" % (k, len(v)) for k, v in o.items()))

    tmp = os.path.join(tmpdir, "gpm_ped.json")
    d9 = build_probe_map(PED_ENTREZ, HG_U133P2, None, tmp, entrez_as_str=True)
    st, dt = commit(tmp, os.path.join(PD, "gene_probe_map_ped.json"))
    report("gene_probe_map_ped.json", st, dt,
           "(%d 基因, 有探针 %d)" % (len(d9), sum(1 for v in d9.values() if v["probes"])))

    if BACKUP_DIR and os.path.isdir(BACKUP_DIR):
        log("\n  被覆盖文件的备份: %s" % os.path.relpath(BACKUP_DIR, ROOT))
    return True


def main():
    ap = argparse.ArgumentParser(
        description="GEO 验证队列原始数据下载 + 8.1/8.2 输入文件构建",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", nargs="?", default="all",
                    choices=["all", "download", "build", "verify"],
                    help="执行阶段（默认 all）")
    ap.add_argument("--force", action="store_true", help="已存在的下载文件也重新下载")
    ap.add_argument("--refresh-entrez", action="store_true",
                    help="用 mygene.info 复核 symbol→Entrez")
    ap.add_argument("--no-backup", action="store_true", help="不备份将被覆盖的文件")
    ap.add_argument("--force-write", action="store_true",
                    help="无条件重写 7 个输入文件（仍先备份）；默认内容一致则不动")
    ap.add_argument("--skip-download", action="store_true", help="等价于 build")
    args = ap.parse_args()

    stage = "build" if (args.skip_download and args.stage == "all") else args.stage

    log("=" * 80)
    log("9.0 验证队列输入准备  |  项目根: %s" % ROOT)
    log("阶段: %s" % stage)
    log("=" * 80)

    setup_backup(not args.no_backup)
    global FORCE_WRITE
    FORCE_WRITE = args.force_write

    # 注释包预检
    for p, tag in [(HG_U219, "hgu219.db (GPL13667)"),
                   (HG_U133P2, "hgu133plus2.db (GPL570)")]:
        log("  %s %-30s %s" % ("✅" if os.path.exists(p) else "❌", tag,
                               os.path.relpath(p, ROOT) if os.path.exists(p) else "缺失"))

    if args.refresh_entrez and stage in ("all", "build"):
        hr("⓪ 用 mygene.info 复核 symbol→Entrez")
        for name, table in [("成人", ADULT_ENTREZ), ("儿童", PED_ENTREZ)]:
            syms = [s for s, e in table.items() if e is not None]
            try:
                got = refresh_entrez(syms)
            except Exception as e:
                log("  %s ✅ 跳过（%s: %s）" % (name, type(e).__name__, e))
                continue
            diff = {s: (table[s], got.get(s)) for s in syms
                    if got.get(s) is not None and got[s] != table[s]}
            log("  %s: 复核 %d 个基因，%s" % (name, len(syms),
                                        "全部一致 ✅" if not diff else "差异 %s" % diff))

    ok = True
    tmpdir = os.path.join(ROOT, "_待清理", "_9.0_tmp")
    os.makedirs(tmpdir, exist_ok=True)

    if stage in ("all", "download"):
        ok = do_download(args) and ok
    if stage in ("all", "build") and ok:
        ok = do_build(args, tmpdir) and ok
    if stage in ("all", "verify") and ok:
        ok = verify() and ok

    # 清理临时目录
    try:
        shutil.rmtree(tmpdir)
    except Exception:
        pass

    hr("衔接：运行 R 端验证")
    log("  cd \"%s\"" % ROOT)
    log("  Rscript 8.1.GSE65682_adult_validation.R       # 期望 n=479 (114 死亡 / 365 存活)")
    log("  Rscript 8.2.GSE66099_pediatric_validation.R   # 期望 181 SepticShock → 108 (13 死亡 / 95 存活)")
    log("  Rscript 9.3.Validation_figures.R              # 出验证图")
    log("")
    log("  结果核对（与正文数值一致即为通过）：")
    log("    成人: 9/9 方向一致 (双侧二项 p=3.91e-03)、多基因评分 p=8.588e-06")
    log("    儿童: 16/20 (双侧 p=1.18e-02)、多基因评分 p=0.03325")
    log("")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
