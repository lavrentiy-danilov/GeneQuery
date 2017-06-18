# python3.5 multi_script.py -i down.txt -o out -t 4 -w 4 -- строка запуска


import os
import time
import shutil
import argparse
import urllib.request
from functools import partial
from multiprocessing import Pool
from subprocess import Popen, PIPE, STDOUT


SPECIES = 'genprime_vM12'
REFERENCE = '/Molly/sakhabeev/reference/{}_kallisto'.format(SPECIES)

def remove(p):
    if os.path.exists(p):
        os.remove(p)

def work_download(url, folder_output):
    name = url.rstrip('/').split('/')[-1]
    path = os.path.join(folder_output, name)
    do_download = True

    if os.path.exists(path):
        size_url = urllib.request.urlopen(url).length
        size_file = os.path.getsize(path)
        if size_url == size_file:
            do_download = False

    if do_download:
        print('[*] Downloading <{}>...'.format(path))
        # def hook(t):
        #     last_b = [0]
        #     def inner(b=1, bsize=1, tsize=None):
        #         if tsize is not None:
        #             t.total = tsize
        #         t.update((b - last_b[0]) * bsize)
        #         last_b[0] = b
        #     return inner
        # with tqdm(unit='B', unit_scale=True, miniters=1, desc=name) as t:
        #     urllib.request.urlretrieve(url, filename=path, reporthook=hook(t))
        cmd_curl = 'curl -L {} -o {}'.format(url, path)
        with Popen(cmd_curl, shell=True, stdout=PIPE, stderr=PIPE) as proc:
            stdout, stderr = proc.communicate()
            if proc.returncode:
                print('[-] Download failed with {} code from {}'.format(proc.returncode, url))
                return None
        print('[+] OK: downloaded <{}>'.format(path))
    else:
        # print('[.] Not downloading <{}>'.format(path))
        pass

    return path


def work_process(path, folder_output, threads):
    tag = os.path.splitext(os.path.basename(path))[0]
    filename_abundance = os.path.join(folder_output, 'abundance', '{}.tsv'.format(tag))

    if os.path.exists(filename_abundance):
        print('[+] Abundance already exists: <{}>'.format(filename_abundance))
    else:
        outdir_fastq_dump = os.path.join(folder_output, 'fastq-dump')
        outdir_fastq_dump_log = os.path.join(folder_output, 'logs', '{}.fastq-dump.log'.format(tag))
        fastq_file_single = os.path.join(outdir_fastq_dump, '{}.fastq'.format(tag))
        fastq_file_first = os.path.join(outdir_fastq_dump, '{}_1.fastq'.format(tag))
        fastq_file_second = os.path.join(outdir_fastq_dump, '{}_2.fastq'.format(tag))

        if os.path.exists(fastq_file_first) and os.path.exists(fastq_file_second):
            print('[.] Using existing paired fastq-dump <{}>'.format(os.path.join(outdir_fastq_dump, '{}_{{1,2}}.fastq'.format(tag))))
            reads = '{} {}'.format(fastq_file_first, fastq_file_second)
            single = ''
        elif os.path.exists(fastq_file_single):
            print('[.] Using existing single fastq-dump <{}>'.format(fastq_file_single))
            reads = fastq_file_single
            single = '--single -l 200 -s 50'
        else:
            print('[*] Fastq-dumping <{}>...'.format(path))
            outdir_fastq_dump_tmp = os.path.join(folder_output, 'fastq-dump_tmp')

            # cmd_fastq_dump = '/Johnny/predeus/programs/sratoolkit.2.5.7-ubuntu64/bin/fastq-dump --split-3 {} --outdir {}'.format(path, outdir_fastq_dump_tmp)
            cmd_fastq_dump = '/Molly/sakhabeev/TOOLS/parallel-fastq-dump-master/parallel-fastq-dump --threads {} --split-3 -s {} --outdir {}'.format(threads, path, outdir_fastq_dump_tmp)

            with open(outdir_fastq_dump_log, 'w') as f_log:
                with Popen(cmd_fastq_dump, shell=True, universal_newlines=True, stdout=f_log, stderr=STDOUT) as proc:
                    proc.wait()
                    if proc.returncode:
                        return 'fastq-dump error {} on {}'.format(proc.returncode, path)

            print('[+] OK: fastq-dump of <{}>'.format(path))

            print('[.] Moving temporary fastq-dump files for <{}>...'.format(tag))
            fastq_file_single_tmp = os.path.join(outdir_fastq_dump_tmp, '{}.fastq'.format(tag))
            fastq_file_first_tmp = os.path.join(outdir_fastq_dump_tmp, '{}_1.fastq'.format(tag))
            fastq_file_second_tmp = os.path.join(outdir_fastq_dump_tmp, '{}_2.fastq'.format(tag))
            if os.path.exists(fastq_file_first_tmp) and os.path.exists(fastq_file_second_tmp):
                shutil.move(fastq_file_first_tmp, fastq_file_first)
                shutil.move(fastq_file_second_tmp, fastq_file_second)
                reads = '{} {}'.format(fastq_file_first, fastq_file_second)
                single = ''
            elif os.path.exists(fastq_file_single_tmp):
                shutil.move(fastq_file_single_tmp, fastq_file_single)
                reads = fastq_file_single
                single = '--single -l 200 -s 50'
            else:
                return 'fastq-dump wierd error - no file after finishing for <{}>'.format(path)

        # -------------

        # print('[*] Removing SRA <{}>...'.format(path))
        # remove(path)

        # -------------

        outdir_kallisto = os.path.join(folder_output, 'kallisto', tag)
        outdir_kallisto_log = os.path.join(folder_output, 'logs', '{}.kallisto.log'.format(tag))
        cmd_kallisto = '/Johnny/predeus/programs/kallisto_linux-v0.43.0/kallisto quant -i {ref} -t {threads} {single} --plaintext -o {outdir} {reads}'.format(ref=REFERENCE, threads=threads, single=single, outdir=outdir_kallisto, reads=reads)

        print('[*] Kallisto <{}>...'.format(tag))
        with open(outdir_kallisto_log, 'w') as f_log:
            with Popen(cmd_kallisto, shell=True, stdout=f_log, stderr=STDOUT) as proc:
                proc.wait()
                if proc.returncode:
                    return 'kallisto error {} on {}'.format(proc.returncode, path)

        print('[+] OK: kallisto on <{}>'.format(tag))

        # -------------

        filename_abundance_from = os.path.join(outdir_kallisto, 'abundance.tsv')
        filename_abundance_to = filename_abundance
        print('[.] Saving adundance to <{}>...'.format(filename_abundance_to))
        # os.rename(filename_abundance_from, filename_abundance_to)
        shutil.copy(filename_abundance_from, filename_abundance_to)

        print('[*] Removing fastq-dump for <{}>...'.format(path))
        remove(os.path.join(outdir_fastq_dump, '{}.fastq'.format(tag)))
        remove(os.path.join(outdir_fastq_dump, '{}_1.fastq'.format(tag)))
        remove(os.path.join(outdir_fastq_dump, '{}_2.fastq'.format(tag)))

        print('[*] Removing kallisto output for <{}>...'.format(tag))
        shutil.rmtree(outdir_kallisto, ignore_errors=True)

        print('[+] OK: {}'.format(path))

    return None


def main():
    _time_start = time.time()

    print('[*] Parsing args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='urls.txt', help='Input filename')
    parser.add_argument('-o', '--output', default='.', help='Output folder')
    parser.add_argument('-w', '--workers', type=int, default=4, help='Workers amount [default: %(default)s]')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Kallisto threads [default: %(default)s]')

    args = parser.parse_args()
    filename_input = args.input
    folder_output = args.output
    workers_amount = args.workers
    threads = args.threads

    os.makedirs(folder_output, exist_ok=True)
    os.makedirs(os.path.join(folder_output, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(folder_output, 'fastq-dump'), exist_ok=True)
    os.makedirs(os.path.join(folder_output, 'fastq-dump_tmp'), exist_ok=True)
    os.makedirs(os.path.join(folder_output, 'kallisto'), exist_ok=True)
    os.makedirs(os.path.join(folder_output, 'abundance'), exist_ok=True)
    if workers_amount == -1 or workers_amount == 0:
        workers_amount = None

    print('[*] Creating pools with <{}> workers each...'.format(workers_amount if workers_amount else 'all'))
    with Pool(workers_amount) as pool_download, Pool(workers_amount) as pool_process, open(filename_input) as f_in:
        try:
            print('[*] Reading urls from <{}>...'.format(filename_input))
            data_urls = (line.rstrip() for line in f_in)

            print('[*] Working...')
            pathes = (path for path in pool_download.imap_unordered(partial(work_download, folder_output=folder_output), data_urls) if path)
            for error in pool_process.imap_unordered(partial(work_process, folder_output=folder_output, threads=threads), pathes):
                if error is not None:
                    print('[!] Error: {}'.format(error))
        except KeyboardInterrupt:
            pool_download.terminate()
            pool_process.terminate()

    print('[+] Done in {:.3f} seconds.'.format(time.time() - _time_start))

if __name__ == '__main__':
main()