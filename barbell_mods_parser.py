# NB: dopo il trimming solo MN, MM e ML vengono ricalcolati (mv/ts/ns/qs rimossi con clean, qs = qualità media pre-trim).
# Gli altri tag di dorado (es. du, e mv/ts/ns/qs se --no_clean) si riferiscono ancora alla read originale non trimmata.
from pysam import FastxFile
from glob import glob
from types import SimpleNamespace
from argparse import ArgumentParser
import os
from concurrent.futures import ProcessPoolExecutor

def unpack_tag(comment : str):
    comment = comment.split('\t')
    res_dict = {}
    types_dict = {}

    for c in comment:
        # TAG:TIPO:VALORE, split massimo 2 volte perché il valore può contenere ':'
        c = c.split(':', 2)
        # salta i campi della descrizione che non sono tag SAM (es. runid=...)
        if len(c) == 3:
            res_dict[c[0]] = c[2]
            types_dict[c[0]] = c[1]

    return SimpleNamespace(**res_dict), types_dict

def parseMM_ML(tag : str):
    # MM: 'A+a.,0,1;C+m?,2;' -> {'A+a.': [0, 1], 'C+m?': [2]}, chiavi nell'ordine del tag
    # ML: 'C,89,108'         -> {'C': [89, 108]}
    tag = tag.rstrip(';').split(';')
    parsed = {}

    for t in tag:
        t = t.split(',')
        parsed[t[0]] = [int(x) for x in t[1:]]

    return parsed

def adjust_tag_trimmed(tag : SimpleNamespace, seq_sine_qua_non : str, bs : int, be : int):
    # seq_sine_qua_non = bp + sequenza trimmata, cioè la read originale fino a be

    mm_tag = parseMM_ML(tag.MM)
    ml_vals = parseMM_ML(tag.ML)['C']

    # posizioni di ogni base sulla read originale
    base_pos = {}
    for pos, base in enumerate(seq_sine_qua_non.upper()):
        base_pos.setdefault(base, []).append(pos)

    new_mm = ''
    new_ml = []
    ml_start = 0

    for key, deltas in mm_tag.items():
        base = key[0]
        codes = key[2:].rstrip('?.')
        # più codici nella stessa sezione (es. C+hm?) -> più valori ML per posizione
        n_codes = 1 if codes.isdigit() else len(codes)
        pos_list = list(range(len(seq_sine_qua_non))) if base == 'N' else base_pos.get(base, [])

        # i delta sono skip: delta -> indice della k-esima occorrenza della base
        idx = []
        cur = -1
        for d in deltas:
            cur += d + 1
            idx.append(cur)

        # ML è un unico array: prendo la fetta di questa sezione
        ml_sec = ml_vals[ml_start:ml_start + len(idx) * n_codes]
        ml_start += len(idx) * n_codes

        # indici oltre la fine di pos_list cadono nella coda tagliata [be, len)
        pos_read_based = [pos_list[i] for i in idx if i < len(pos_list)]
        to_keep = [k for k, p in enumerate(pos_read_based) if (p >= bs) and (p < be)]

        # indici locali al frammento: tolgo le occorrenze nel prefisso
        n_pre = sum(1 for p in pos_list if p < bs)
        local = [idx[k] - n_pre for k in to_keep]
        new_deltas = [l - prev - 1 for prev, l in zip([-1] + local[:-1], local)]

        new_mm += key + ''.join(f',{d}' for d in new_deltas) + ';'
        for k in to_keep:
            new_ml.extend(ml_sec[k * n_codes:(k + 1) * n_codes])

    assert ml_start == len(ml_vals), 'ML non coerente con MM'

    return new_mm, 'C' + ''.join(f',{v}' for v in new_ml)

def write_fastq_record(fh, name, seq, qual, tag, types):
    fields = [f"{k}:{types[k]}:{v}" for k, v in vars(tag).items()]
    header = '\t'.join([name] + fields)      # tab, non spazi: serve per minimap2 -y
    fh.write(f"@{header}\n{seq}\n+\n{qual}\n")

def process_file(f, outdir, clean):
    out_path = os.path.join(outdir, os.path.basename(f))

    with FastxFile(f) as fh, open(out_path, 'w') as out:
        for entry in fh:
            tag, types_dict = unpack_tag(entry.comment)

            assert getattr(tag, 'bf', '0') == '0'
            bs, be = int(tag.bs), int(tag.be)
            prefix = getattr(tag, 'bp', '')   # assente = niente tagliato all'inizio
            assert len(prefix) == bs and bs + len(entry.sequence) == be

            #handle possibly mismatching tags after trimming
            tag.MN = len(entry.sequence)
            types_dict['MN'] = 'i'
            if clean:
                for key in ('mv', 'ts', 'ns', 'qs'):
                    vars(tag).pop(key, None)
                    types_dict.pop(key, None)


            new_mm, new_ml = adjust_tag_trimmed(tag, prefix + entry.sequence, bs, be)
            tag.ML = new_ml
            tag.MM = new_mm
            # bp manca quando bs = 0
            for key in ('bs', 'be', 'bf', 'bp'):
                vars(tag).pop(key, None)

            write_fastq_record(out, entry.name, entry.sequence, entry.quality, tag, types_dict)
    print('done file: ' + f)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--fastqs_folder", required=True)
    parser.add_argument("--out_folder", required=True)
    parser.add_argument("--cores", required=False, default=None, type=int)
    # pulizia di mv/ts/ns/qs attiva di default, --no_clean per disattivarla
    parser.add_argument("--no_clean", dest='clean', action='store_false')
    args = parser.parse_args()

    fastqs_folder = args.fastqs_folder
    outdir = args.out_folder
    cores = args.cores
    clean = args.clean

    fqs = glob(fastqs_folder + '/*.fastq')

    if cores is None:
        for f in fqs:
            process_file(f, outdir, clean)
    else:
        with ProcessPoolExecutor(max_workers=cores) as executor:
            # list() serve a far risalire le eccezioni dei worker
            list(executor.map(process_file, fqs, [outdir] * len(fqs), [clean] * len(fqs)))
