#!/usr/bin/env python3

from math import log2
import re
import sys
from pathlib import Path

BENCH_RE = re.compile(r"^Benchmarking\s+(.+?)/T(\d+)_H(\d+)")
TIME_RE = re.compile(r"time:\s*\[([^\]]+)\]")
COLLECT_RE = re.compile(r"Collecting\s+(\d+)\s+samples.*\(([^)]+)\s+iterations\)")
PROOF_RE = re.compile(r"proof.*?size.*?([0-9]+(?:[.,][0-9]+)?)", re.IGNORECASE)

def main():
    if len(sys.argv) < 2:
        print("Usage: gen_plots.py <input-criterion-bench-file>")
        sys.exit(1)
    input_path = Path(sys.argv[1])
    
    text = input_path.read_text(encoding='utf-8')
    lines = text.splitlines()

    results = []

    # Use a mapping keyed by (group,bench_id,T,H) to accumulate info because
    # the "Benchmarking ..." header appears multiple times (warming, collecting,
    # analyzing) and we want to aggregate fields that may appear on different
    # lines.
    in_progress = {}
    last_key = None
    pending_proof_bytes = None

    for line in lines:
        line = line.rstrip('\n')
        m = BENCH_RE.match(line)
        if m:
            bench_path, t_str, h_str = m.groups()
            parts = bench_path.split('/')
            group = parts[0]
            key = (group, int(t_str), int(h_str))
            last_key = key
            if key not in in_progress:
                in_progress[key] = {
                    'bench_group': group,
                    'T': int(t_str),
                    'H': int(h_str),
                    'median': None,
                    'unit': None,
                    'proof_size_bytes': None,
                }
        if last_key is None:
            continue
        mtime = TIME_RE.search(line)
        if mtime:
            bracket = mtime.group(1).strip()
            parts = bracket.split()
            median =  parts[0]
            if parts[1] == 'µs':
                unit = 'us'
            else:
                unit = parts[1]
            in_progress[last_key]['median'] = median
            in_progress[last_key]['unit'] = unit
            continue
        mcol = COLLECT_RE.search(line)
        if mcol:
            # finish block only if time info already present
            if in_progress[last_key].get('median') is not None:
                results.append(in_progress.pop(last_key))
            continue
        # try match proof size
        mproof = PROOF_RE.search(line)
        if mproof:
            num = mproof.group(1)
            # normalize number
            num = num.replace(',', '.')
            try:
                val = float(num)
            except Exception:
                val = None
            if val is not None:
                # if there is an active in_progress entry for last_key, attach immediately
                if last_key in in_progress:
                    in_progress[last_key]['proof_size_bytes'] = val
            continue
    # flush any remaining completed entries
    for key, val in list(in_progress.items()):
        if val.get('median') is not None:
            results.append(val)

    by_t_h = {}
    for r in results:
        key = (r['T'], r['H'])
        by_t_h.setdefault(key, []).append(r)

    # helper: convert median string + unit -> milliseconds (float)
    def median_to_ms(median_str, unit):
        try:
            val = float(median_str)
        except Exception:
            return None
        u = (unit or '').lower()
        if u == 'us' or u == '\u00b5s':
            return val * 1e-3
        if u == 'ns':
            return val * 1e-6
        if u == 'ms':
            return val
        if u == 's':
            return val * 1e3
        # unknown unit: assume seconds
        return val * 1e3

    # Build a mapping per (T,H) of bench_group -> median_ms
    per_param = {}
    for (t, h), items in by_t_h.items():
        grp = {}
        proof_bytes = None
        for it in items:
            g = it.get('bench_group')
            med = it.get('median') or it.get('median_s') or it.get('median_ms')
            unit = it.get('unit')
            ms = median_to_ms(med, unit)
            grp[g] = ms
            if proof_bytes is None and it.get('proof_size_bytes') is not None:
                proof_bytes = it.get('proof_size_bytes')
        if proof_bytes is not None:
            grp['proof_size_bytes'] = proof_bytes
        per_param[(t, h)] = grp

    # Create output directory for per-T datapoints
    out_dir = input_path.parent / 'bench_data'
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {}
    # For each T, create files with rows: time_ms H for baseline, prover, verifier
    all_Ts = sorted({k[0] for k in per_param.keys()})
    for t in all_Ts:
        baseline_points = []
        prover_points = []
        compiled_verifier_points = []

        for (tt, h), grp in per_param.items():
            h = log2(h)
            if tt != t:
                continue
            # raw components
            ve = grp.get('vector_encryption')
            vc = grp.get('verifier_computation')
            rd = grp.get('result_decryption')
            pm = grp.get('plaintext_matvec')
            pr = grp.get('prover_computation')

            # baseline: plaintext_matvec
            if pm is not None:
                baseline_points.append((pm, h))

            # prover computation
            if pr is not None:
                prover_points.append((pr, h))

            if ve is not None and vc is not None and rd is not None:
                compiled_verifier_points.append((ve + vc + rd, h))

        # helper to write a list of points to a file
        def write_points(name, pts):
            pts_sorted = sorted(pts, key=lambda x: x[1])
            out_f = out_dir / f"T{t}__{name}.dat"
            with out_f.open('w') as fh:
                fh.write('# time_ms H\n')
                for time_ms, h in pts_sorted:
                    fh.write(f"{time_ms:.6f} {h}\n")
            return pts_sorted
        
        summary[t] = {
            'baseline': write_points('baseline', baseline_points),
            'compiled_verifier': write_points('compiled_verifier', compiled_verifier_points),
        }

    # Generate per-T PGFPlots figure fragments (.tex), one file per T.
    figs_dir = input_path.parent / 'bench_figs'
    figs_dir.mkdir(parents=True, exist_ok=True)

    # Generate LaTeX table for prover times: columns = T values, rows = lg n values
    tables_dir = input_path.parent / 'bench_tables'
    tables_dir.mkdir(parents=True, exist_ok=True)

    # gather all unique H values and sort by log2(H)
    all_H = sorted({k[1] for k in per_param.keys()}, key=lambda x: log2(x))
    all_T = sorted({k[0] for k in per_param.keys()})

    table_path = tables_dir / 'prover_times.tex'
    with table_path.open('w') as tf:
        tf.write('% Auto-generated LaTeX table of prover times (seconds)\n')
        cols = '|c|' + 'c' * len(all_T) + '|'
        tf.write('\\begin{table}\n')
        tf.write('\\centering\n')
        tf.write('\\begin{tabular}{%s}\n' % cols)
        tf.write('\\hline\n')
        header_cells = ['\\diagbox{$\\lg n$}{$\\lg t$}'] + [f"${int(log2(t))}$" for t in all_T]
        tf.write(' & '.join(header_cells) + ' \\\\ \n\\hline\n')

        for H in all_H:
            lg = log2(H)
            lg_label = f"{int(lg)}" if abs(lg - round(lg)) < 1e-9 else f"{lg:.2f}"
            row = [lg_label]
            for t in all_T:
                grp = per_param.get((t, H), {})
                pr_ms = grp.get('prover_computation')
                if pr_ms is None:
                    cell = '-' 
                else:
                    pr_s = pr_ms / 1000.0
                    cell = f"{pr_s:.2f}"
                row.append(cell)
            tf.write(' & '.join(row) + ' \\\\\n')
        tf.write('\\hline\n')

        tf.write('\\end{tabular}\n')
        tf.write('\\caption{Prover times (in seconds).}\n')
        tf.write('\\end{table}\n')

    print(f"Wrote prover times LaTeX table to {table_path}")

    proof_table_path = tables_dir / 'proof_sizes.tex'
    with proof_table_path.open('w') as pf:
        pf.write('% Auto-generated LaTeX table of proof sizes (kB)\n')
        cols = '|c|' + 'c' * len(all_T) + '|'
        pf.write('\\begin{table}\n')
        pf.write('\\centering\n')
        pf.write('\\begin{tabular}{%s}\n' % cols)
        pf.write('\\hline\n')
        header_cells = ['\\diagbox{$\\lg n$}{$\\lg t$}'] + [f"${int(log2(t))}$" for t in all_T]
        pf.write(' & '.join(header_cells) + ' \\\\ \\hline\n')

        for H in all_H:
            lg_label = str(int(log2(H)))
            row = [lg_label]
            for t in all_T:
                grp = per_param.get((t, H), {})
                ps_bytes = grp.get('proof_size_bytes')
                ps_kb = float.__round__(ps_bytes / 1000)
                cell = ps_kb
                row.append(str(cell))
            pf.write(' & '.join(row) + ' \\\\\n')

        pf.write('\\hline\n')
        pf.write('\\end{tabular}\n')
        pf.write('\\caption{Proof sizes (in kB).}\n')
        pf.write('\\end{table}\n')

    print(f"Wrote proof sizes LaTeX table to {proof_table_path}")

    for t, _ in summary.items():
        tex_path = figs_dir / f"fig_T{t}_verifier.tex"
        with tex_path.open('w') as fh:
            fh.write('% Auto-generated PGFPlots figure for T={0}\n'.format(t))
            fh.write('\\begin{figure}[ht]\n')
            fh.write('\\centering\n')
            fh.write('\\begin{tikzpicture}\n')
            fh.write('  \\begin{axis}[\n')
            fh.write('    xlabel={$\\lg n$},\n')
            fh.write('    ylabel={Time (ms)},\n')
            fh.write('    grid=major,\n')
            fh.write('    legend style={at={(0,1)}, anchor=north west, legend cell align=left}\n')
            fh.write('  ]\n')

            # For each metric, if the corresponding file exists and has data, add a plot line
            metric_names = ['baseline', 'compiled_verifier']
            labels = {'baseline': 'Baseline', 'compiled_verifier': 'Compiled verifier'}
            colors = {'baseline': 'blue', 'compiled_verifier': 'orange'}

            for m in metric_names:
                data_file = out_dir / f"T{t}__{m}.dat"
                if data_file.exists():
                    # use table reading; note column order: time_ms H
                    fh.write('    \\addplot+[mark=*, thick, color={0}] table[x index=1,y index=0] {{{1}}};\n'.format(colors[m], str(data_file)) )
                    fh.write('    \\addlegendentry{{{0}}}\n'.format(labels[m]))

            fh.write('  \\end{axis}\n')
            fh.write('\\end{tikzpicture}\n')
            fh.write('\\caption{Timing for $t = 2^{%s}$ across $n$ values.}\n' % int(log2(t)))
            fh.write('\\end{figure}\n')

    print(f"Wrote PGFPlots fragments to {figs_dir}")
if __name__ == '__main__':
    main()
