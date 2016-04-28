#binaries
#organism can be ['human', 'chicken', 'yeast']
organism = "human"

#assembly can be ['trinity', 'oases']
assembly = "trinity"

#############################################
# remember sailfish binary is from bin NOT from build/src
#sailfish_binary = "/home/hirak/RapCompress/sailfish/bin/sailfish"
sailfish_binary = "/home/hirak/quark/sailfish/build/src/sailfish"
############################################
sailfish_threads = 20
data_path = "/mnt/scratch1/hirak/RapCompressData/"
index_path = data_path + "sailfish/sailfish_index/index_31/"
ref_path = data_path + "flux/transcriptome/ref.transcripts.fa"

rule make_sailfish_index:
    input:
        ref = ref_path
    output:
        index = index_path
    run:
        shell("mkdir -p {}".format(output.index))
        shell("{} index -p {} -t {} -o {} -k 31".format(sailfish_binary, sailfish_threads,input.ref, output.index))

rule run_sailfish:
    input:
        index = data_path + "sailfish/sailfish_index/index_31/"
    output:
        sfpath = data_path + "quarksailfish/sailfish_quant"
    threads:
        sailfish_threads
    run:
        shell ("mkdir -p {output.sfpath}")
        one = data_path + "flux/reads/r1.fq.gz"
        two = data_path + "flux/reads/r2.fq.gz"
        shell("{} quant -i {} -l IU -p {} -1 <(gunzip -c {}) -2 <(gunzip -c {}) -o {} --dumpEq".format(sailfish_binary, input.index, sailfish_threads, one, two, output.sfpath))


