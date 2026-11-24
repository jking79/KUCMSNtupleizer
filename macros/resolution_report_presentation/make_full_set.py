#!/usr/bin/env python3
"""Build the full plot-family gallery from report-selected resolution fits."""
import collections
import html
import re

import make_prototypes as style  # Also refreshes the three pinned review prototypes.


FAMILIES = [
    ("family_01_uncal_gain1", "Uncalibrated · Gain1Only", "SRO, DRO, and ZEE, split by EB and EE; uncalibrated and Gain1Only are the defaults."),
    ("family_02_gain1_chronus_vs_default", "Gain1Only · Chronus versus default", "SRO and DRO Chronus-calibrated and uncalibrated comparisons, split by EB and EE."),
    ("family_03_allgains_chronus_vs_default", "AllGains · Chronus versus default", "SRO and DRO Chronus-calibrated and uncalibrated comparisons using AllGains, split by EB and EE."),
    ("family_04_chronus_gain1_vs_allgains", "Chronus · Gain1Only versus AllGains", "SRO and DRO Chronus-calibrated comparisons, split by EB and EE."),
    ("family_05_chronus_allgains_eb_ee", "Chronus · AllGains · EB versus EE", "SRO and DRO comparisons across EB and EE."),
    ("family_06_ecal_dpg_eb", "ECAL DPG · EB SRO/DRO/ZEE", "EB-only, uncalibrated, Gain1Only SRO, DRO, and ZEE comparison."),
]


def era_description(tag):
    match = re.search(r"(?:eg_)(16|17|18|22|23|24|25|26)", tag)
    if not match:
        raise ValueError(f"Could not determine year from dataset tag: {tag}")
    short = match.group(1)
    year = 2000 + int(short)
    if year <= 2018:
        era = f"{year} UL"
        if "_ABC_" in tag:
            era += " ABC"
        elif "_D_" in tag:
            era += " D"
        production = "MiniAODv2"
    elif year in (2022, 2023):
        era = str(year)
        production = "19Dec2023"
    elif year == 2024:
        era = str(year)
        production = "MINIv6NANOv15"
    else:
        era = f"{year} Prompt"
        production = "Prompt"
    first_line = era if production == "Prompt" else f"{era} · {production}"
    slug = re.sub(r"[^A-Za-z0-9]+", "_", era).strip("_")
    return year, era, production, first_line, slug


def load_eras():
    eras = collections.defaultdict(dict)
    for hist, fitrow in style.audit.items():
        if fitrow.get("processing_status") != "success":
            continue
        era_key = (int(fitrow["run_start"]), int(fitrow["run_end"]), fitrow["dataset_tag"])
        bin_files = list(style.REPORT.glob(f"{hist}_NSC_resfit_resolution_bin_results.csv"))
        if len(bin_files) != 1:
            raise RuntimeError(f"Expected one bin-results CSV for {hist}, found {len(bin_files)}")
        key = (fitrow["pair_type"], fitrow["gain_selection"], fitrow["calibration_state"], fitrow["region"])
        eras[era_key][key] = (fitrow, style.read_rows(bin_files[0]))
    return eras


def label_pair_calibration(pair, cali):
    return pair + (" Chronus" if cali == "Cali" else "")


def family_plots(era_records):
    plots = []

    def add(slug, region, subtitle, specs, xmax):
        curves = []
        for key, label, pair, gain, cali, reg in specs:
            record = era_records.get((pair, gain, cali, reg))
            if record is None:
                return
            curves.append((key, label, record))
        plots.append((slug, region, subtitle, curves, xmax))

    for region in ("EB", "EE"):
        add("family_01_uncal_gain1", region, region,
            [(pair, pair, pair, "Gain1Only", "NoCali", region) for pair in ("SRO", "DRO", "ZEE")], 1600)
        add("family_02_gain1_chronus_vs_default", region, region,
            [(f"{pair}_{cali}", label_pair_calibration(pair, cali), pair, "Gain1Only", cali, region)
             for pair in ("SRO", "DRO") for cali in ("Cali", "NoCali")], 1600)
        add("family_03_allgains_chronus_vs_default", region, region,
            [(f"{pair}_{cali}", label_pair_calibration(pair, cali), pair, "AllGains", cali, region)
             for pair in ("SRO", "DRO") for cali in ("Cali", "NoCali")], 3200)
        add("family_04_chronus_gain1_vs_allgains", region, region,
            [(f"{pair}_{gain}", pair + (" AllGains" if gain == "AllGains" else ""), pair, gain, "Cali", region)
             for pair in ("SRO", "DRO") for gain in ("Gain1Only", "AllGains")], 3200)
    add("family_05_chronus_allgains_eb_ee", "EB_EE", "Chronus · AllGains",
        [(f"{pair}_{region}_Cali", f"{pair} {region}", pair, "AllGains", "Cali", region)
         for pair in ("SRO", "DRO") for region in ("EB", "EE")], 3200)
    add("family_06_ecal_dpg_eb", "EB", "EB",
        [(pair, pair, pair, "Gain1Only", "NoCali", "EB") for pair in ("SRO", "DRO", "ZEE")], 800)
    return plots


def html_caption(era_label, subtitle, curves):
    parts = [f"{era_label} · {subtitle}" if subtitle != era_label else era_label]
    excluded = []
    for _, label, record in curves:
        count = int(record[0].get("excluded_tail_bins", 0))
        if count:
            excluded.append(f"{label} ({count})")
    if excluded:
        parts.append("High-energy points inconsistent with the fit body were excluded for " + ", ".join(excluded) + "; they remain visible as same-color crosses.")
    return "<br>".join(html.escape(part) for part in parts)


def build_gallery(generated, eras):
    grouped = collections.defaultdict(lambda: collections.defaultdict(list))
    for item in generated:
        grouped[item["family"]][item["era_key"]].append(item)
    nav = ['<nav class="family-nav"><strong>Full set · jump to family:</strong>']
    sections = ['<!-- FULL_SET_GALLERY_START -->', '<h2>Full plot set · grouped by family</h2>',
                '<p>Each family contains plots for all complete Run 2 and Run 3 era records in the fit report. Production labels appear inside the plot. Tail-excluded high-energy points retain their series colors and use crosses; captions identify the fits where this occurred.</p>']
    for slug, title, description in FAMILIES:
        nav.append(f'<a href="#{slug}">{html.escape(title)}</a>')
        family_items = grouped.get(slug, {})
        total = sum(len(v) for v in family_items.values())
        sections.append(f'<section class="plot-family" id="{slug}"><h3>{html.escape(title)} <span class="tag">{total} plots</span></h3><p>{html.escape(description)}</p>')
        for era_key in sorted(family_items, key=lambda k: (int(k[0]), k[2])):
            year, era_label, production, first_line, era_slug = era_description(era_key[2])
            items = family_items[era_key]
            sections.append(f'<details class="era-group"><summary>{html.escape(first_line)} · {len(items)} plot(s)</summary><div class="grid">')
            for item in items:
                rel = item["relative_path"]
                name = html.escape(item["title"])
                caption = html_caption(item["era_label"], item["subtitle"], item["curves"])
                sections.append(f'<figure><a href="{rel}.pdf"><img loading="lazy" src="{rel}.png" alt="{name}"></a><figcaption><strong>{name}</strong><br>{caption}<br><a href="{rel}.pdf">PDF</a></figcaption></figure>')
            sections.append('</div></details>')
        sections.append('</section>')
    nav.append('</nav>')
    sections.append('<!-- FULL_SET_GALLERY_END -->')
    page = style.HERE / "index.html"
    text = page.read_text()
    start, end = '<!-- FULL_SET_GALLERY_START -->', '<!-- FULL_SET_GALLERY_END -->'
    if start in text and end in text:
        before = text.split(start)[0]
        after = text.split(end, 1)[1]
        text = before + after
    marker = '<h2>Working plots · reference only</h2>'
    css = '.family-nav{position:sticky;top:0;z-index:5;display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:10px 12px;background:#fff;border:1px solid var(--line);border-radius:8px;margin:20px 0}.family-nav a{color:var(--blue);font-size:14px}.plot-family{padding-top:8px}.plot-family h3{font-size:19px;margin:26px 0 5px}.era-group{margin:10px 0;border:1px solid var(--line);border-radius:8px;background:#fff}.era-group summary{cursor:pointer;padding:11px 14px;font-weight:650}.era-group .grid{padding:0 12px 12px}.plot-family .tag{float:right}.plot-family p{margin-top:4px}'
    text = text.replace('</style>', css + '\n</style>')
    block = '\n'.join(sections[:3] + nav + sections[3:]) + '\n'
    if marker not in text:
        raise RuntimeError("Could not find gallery insertion point in index.html")
    text = text.replace(marker, block + marker, 1)
    page.write_text(text)


def main():
    eras = load_eras()
    generated = []
    skipped = []
    for era_key, era_records in sorted(eras.items(), key=lambda item: (item[0][0], item[0][2])):
        year, era_label, production, first_line, era_slug = era_description(era_key[2])
        run_text = "Run 2 (13 TeV)" if era_key[0] < 350000 else "Run 3 (13.6 TeV)"
        for family_slug, family_title, _ in FAMILIES:
            for plot_slug, region, second_line, curves, xmax in family_plots(era_records):
                if plot_slug != family_slug:
                    continue
                filename = f"{family_slug}/{era_slug}_{region}"
                subtitle = f"{first_line}\n{second_line}"
                style.make_plot(filename, curves, subtitle, run_text, xmax_override=xmax)
                generated.append({
                    "family": family_slug,
                    "era_key": era_key,
                    "era_label": era_label,
                    "subtitle": second_line,
                    "title": f"{era_label} · {second_line}",
                    "relative_path": f"plots/{filename}",
                    "curves": curves,
                })
    build_gallery(generated, eras)
    print(f"Generated {len(generated)} full-set plots across {len(FAMILIES)} families and {len(eras)} era records.")


if __name__ == "__main__":
    main()
