import pod5 as p5
import matplotlib.pyplot as plt
import numpy as np
import pysam


# Pfad zur POD5-Datei
POD5_FILE = "test_files/example_read.pod5"
BAM_FILE = "test_files/filtered.bam"
MID_BASES = 20


# === Events mit Sequenz aus BAM extrahieren ===
def get_all_segments_with_seq(bam_path, read_id, include_seq=False):
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for aln in bam.fetch(until_eof=True):
            if str(read_id) in aln.query_name and aln.has_tag("mv"):
                seq = aln.query_sequence
                moves = list(aln.get_tag("mv"))
                segments = []
                sample_index = 0
                base_index = 0
                for m in moves:
                    sample_index += m
                    if m > 0 and base_index < len(seq):
                        if include_seq:
                            segments.append((sample_index - m, sample_index, base_index, seq))
                        else:
                            segments.append((sample_index - m, sample_index, seq[base_index]))
                        base_index += 1
                return segments
    return []

# === Signal aus POD5 lesen ===
with p5.Reader(POD5_FILE) as reader:
    reads_iter = reader.reads()
    for i, read in enumerate(reads_iter):
        if i == 9:
            read_id = str(read.read_id)
            signal = read.signal_pa
            rate = read.run_info.sample_rate
            time = np.arange(len(signal)) / rate
            break

# === Event-Segmente ermitteln
segments_raw = get_all_segments_with_seq(BAM_FILE, read_id, include_seq=True)

# Schneide auf Mitte
n = len(segments_raw)
half = MID_BASES // 2
start_idx = max(0, n // 2 - half)
end_idx = start_idx + MID_BASES
cut_segments = segments_raw[start_idx:end_idx]

# Dann erst gültige 5-mer extrahieren
plot_segments = []
for s, e, base_index, seq in cut_segments:
    start = base_index - 2
    end = base_index + 3
    if start >= 0 and end <= len(seq):
        five_mer = seq[start:end]
        plot_segments.append((s, e, five_mer))

if not plot_segments:
    raise RuntimeError("Keine 5-mer-Events in Auswahl")


# === Signal-Ausschnitt laden
start_sample = plot_segments[0][0]
end_sample = plot_segments[-1][1]
signal = signal[start_sample:end_sample]
time = time[start_sample:end_sample]

# === Plot erstellen
plt.figure(figsize=(14, 5))

# Basenfarben nach der **ersten Base** des 5-mer
base_colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}

for s, e, five_mer in plot_segments:
    s_rel = s - start_sample
    e_rel = e - start_sample
    if s_rel < 0 or e_rel > len(signal):
        continue

    t_vals = time[s_rel:e_rel]
    y_vals = signal[s_rel:e_rel]
    color = base_colors.get(five_mer[0], 'gray')

    # Streupunkte einfärben
    plt.scatter(t_vals, y_vals, color=color, s=20, alpha=0.7, zorder=3)

    # Label direkt knapp oberhalb des letzten Punkts platzieren
    x_center = (t_vals[0] + t_vals[-1]) / 2
    y_label = np.mean(y_vals) + 0.5

    first, rest = five_mer[0], five_mer[1:]

    # Zeichne zuerst den eingefärbten Teil
    plt.text(x_center, y_label, first,
             ha='right', fontsize=12, zorder=4,
             bbox=dict(facecolor=color, alpha=0.3, edgecolor=color, boxstyle='round,pad=0.15'))

    # Dann den Rest direkt daneben (schwarzer Text, kein Rahmen)
    plt.text(x_center + 0.000001, y_label, rest,  # kleiner X-Offset
             ha='left', fontsize=12, color='black', zorder=4)

plt.xlabel("Time (s)")
plt.ylabel("Signal (pA)")
plt.title(f"{MID_BASES} Events from read {read_id}")
plt.tight_layout()
plt.subplots_adjust(top=0.87)
plt.show()
