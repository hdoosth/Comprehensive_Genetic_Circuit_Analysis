# Load packages
import csv
from pprint import pprint as pp
import json
import os
import pandas as pd
import numpy as np
# import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import datetime
import json
# import glob
from matplotlib import gridspec
import shlex, subprocess

default_prefix = './pavana/input/'
default_output = './pavana/output/'
# with open(".\pavana\input//15B_R_1.fastq",'rU') as ins:
#     count = 0
#     x=[]
#     for lines in ins:
#         count +=1
#         field = lines.strip().split()
#         # if lines[0] == '@':
#         #     if field[0] == '@SQ':
#         #         print(field)
#         # else:
#             # x.append(field[1])
#         print(lines)
#         if count > 100:
#             break
#     print(set(x),len(x))
#     for i in set(x):
#         print(i, x.count(i))

def create_profile(file_name, prefix=default_prefix, output_path=default_output, sample_id=""):

    if sample_id =="":
        sample_id = file_name.split('.')[0]
    SAM_file = prefix + file_name
    circuit_fwd_avg_read = {}
    circuit_rev_avg_read = {}
    reads = {}
    profile = {}
    genomes = {}
    flag_genomecreation = 0
    nline = 0
    mappedline = 0
    min_length = 10
    stats = {}
    max_values = {}
    with open(SAM_file, 'rU') as ins:
        for lines in ins:
            field = lines.strip().split()
            if lines[0] == '@':
                if field[0] == '@SQ':
                    name = field[1].split(':')
                    length = field[2].split(':')
                    genomes[name[1]] = length[1]
                    profile[name[1] + '_fwd'] = np.zeros(int(length[1]))
                    profile[name[1] + '_rev'] = np.zeros(int(length[1]))
                    reads[name[1] + '_fwd'] = {}
                    reads[name[1] + '_rev'] = {}
            if lines[0] != '@':
                if field[2] in genomes and int(field[4]) >= 10:                     
                    start = min(int(field[3]),int(field[7]))
                    readlen = abs(int(field[8]))-1
                    if int(field[1]) == 83 or int(field[1]) == 163:
                        if readlen >= min_length:
                            if start not in reads[field[2]+'_fwd']:
                                reads[field[2]+'_fwd'][start] = []
                            reads[field[2]+'_fwd'][start].append(readlen)
                            profile[field[2]+'_fwd'][start:start+readlen] += np.ones(readlen)
                            mappedline += 1
                            print(lines)
                    elif int(field[1]) == 99 or int(field[1]) == 147:
                        if readlen >= min_length:
                            if start not in reads[field[2]+'_rev']:
                                reads[field[2]+'_rev'][start] = []
                            reads[field[2]+'_rev'][start].append(readlen)
                            profile[field[2]+'_rev'][start:start+readlen] += np.ones(readlen)
                            mappedline += 1

                nline += 1
    for seq, length in genomes.items():
        max_values[seq+'_fwd'] = max(profile[seq+'_fwd'])
        max_values[seq+'_rev'] = max(profile[seq+'_rev'])
    print(str(nline) + " reads found in SAM file " + str(mappedline) + " met qc to be mapped into profile")
    stats['reads_sam'] = nline
    stats['reads_profile'] = mappedline
    stats['max_values'] = max_values
    stats['genomes'] = genomes
    print("Profile Generation Complete - Saving")
    json_reads = json.dumps(reads)
    with open(output_path+sample_id+"_reads.json","w") as f_reads:
        f_reads.write(json_reads)
    
    profile_for_save = {}
    for key, length in genomes.items():
        profile_for_save[key + '_fwd'] = profile[key + '_fwd'].tolist()
        profile_for_save[key + '_rev'] = profile[key + '_rev'].tolist()
    json_profile = json.dumps(profile_for_save)
    with open(output_path+sample_id+"_profile.json","w") as f_profile:
        f_profile.write(json_profile)
    print("Save Complete")
    
    return profile, stats

def load_profile(file_name, genomes, output_path=default_output):
    profile = {}
    for key, length in genomes.items():
        profile[key + '_fwd'] = np.zeros(length)
        profile[key + '_rev'] = np.zeros(length)
    with open(output_path+file_name,"rU") as f_profile:
        profile_load = json.load(f_profile)
    for key, length in genomes.items():
        profile[key + '_fwd'] = np.array(profile_load[key + '_fwd'])
        profile[key + '_rev'] = np.array(profile_load[key + '_rev'])
    return profile

def load_trna_dict(filename):
	gene_dict = {}
	lines = [open(filename, 'r').read().strip("\n")][0].split('\n')
	for line in lines:
		if not line.startswith('#'):
			tokens = line.split('\t')
			genome, feature, start, end, strand = tokens[0], tokens[2], int(tokens[3]), int(tokens[4]), tokens[6]
			if feature in ['tRNA','rRNA']:
				name = (tokens[-1].split('Name='))[-1].split(';')[0]
				if genome not in gene_dict:
					gene_dict[genome] = {}
				if name not in gene_dict[genome]:
					gene_dict[genome][name] = {'start': start, 'end': end, 'strand': strand, 'feature': feature}
				if name in gene_dict[genome]:
					gene_dict[genome][name+'_v'] = {'start': start, 'end': end, 'strand': strand, 'feature': feature}
	return gene_dict

def normalize(profile_raw, sample_id, gff_file, prefix=default_prefix, output_path=default_output):
    trRNA_dict = load_trna_dict(prefix + gff_file)
    profile = {} # normalized profile values
    profile_for_save = {}
    for key in profile_raw.keys():
        total_reads = 0
        if key[-4::] == '_fwd':
            genome = key[0:-4]
            total_reads += sum(profile_raw[genome + '_fwd']) + sum(profile_raw[genome + '_rev'])
            if genome in trRNA_dict.keys():
                for entry in trRNA_dict[genome]:
                    if not entry.endswith('_v'):
                        # print(trRNA_dict[genome][entry])
                        strand_name = '_fwd' if trRNA_dict[genome][entry]['strand'] =='+' else '_rev'
                        total_reads -= sum(profile_raw[genome+strand_name][trRNA_dict[genome][entry]['start']:trRNA_dict[genome][entry]['end']])
    for key in profile_raw.keys():
        profile[key] = profile_raw[key]/total_reads*1e9
        profile_for_save[key] = profile[key].tolist()
    json_profile = json.dumps(profile_for_save)
    with open(output_path+sample_id+"_profile_norm.json","w") as f_profile_norm:
        f_profile_norm.write(json_profile)
    print('Normalized profile created and saved for:', sample_id)
    return profile

def plot_regions_of_interest(profile, regions_of_int, output_path=default_output):
    for chrom in regions_of_int.keys():
        for cur_region in regions_of_int[chrom]:
                        
            fig = plt.figure(figsize=(12.0,4.0))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.2])
            ax = plt.subplot(gs[0])

            avg_profile_fwd = np.zeros(len(range(cur_region[0],cur_region[1])))
            avg_profile_rev = np.zeros(len(range(cur_region[0],cur_region[1])))
            for key in profile.keys():
                if key[-4::] == '_fwd':
                    genome = key[0:-4]
                    print(genome,chrom)
                    if genome == chrom:
                        super_threshold_indices = profile[genome + '_fwd'] < 1
                        profile[genome + '_fwd'][super_threshold_indices] = 1
                        super_threshold_indices = profile[genome + '_rev'] < 1
                        profile[genome + '_rev'][super_threshold_indices] = 1
                        # print(len(avg_profile_fwd))
                        # print(len(profile[genome + '_fwd'][cur_region[0]:cur_region[1]]))
                        avg_profile_fwd += np.log10(profile[genome + '_fwd'][cur_region[0]:cur_region[1]])
                        avg_profile_rev += -np.log10(profile[genome + '_rev'][cur_region[0]:cur_region[1]])
            ax.fill_between(np.array(range(cur_region[0],cur_region[1])), 0.0, avg_profile_fwd, facecolor='grey', linewidth=0.0)
            ax.fill_between(np.array(range(cur_region[0],cur_region[1])), 0.0, avg_profile_rev, facecolor='lightpink', linewidth=0.0)
            ax.set_xlim(cur_region)
            ax.ticklabel_format(style='plain')
            ax.set_ylim([-5,5])
            for axis in ['bottom','top','right']:
                ax.spines[axis].set_linewidth(1.0)
                ax.spines[axis].set_color('k')
            ax.spines['left'].set_linewidth(1.0)
            ax.spines['left'].set_color('k')
                        
            ax.tick_params(axis='y', which='major', labelsize=7, length=4, width=0.25, direction='out')
            ax.tick_params(axis='x', which='major', labelsize=7, length=4, width=0.25, direction='out')
            ax.set_yticks([-5,0,5])
            ax.set_yticklabels([-5,0,5])

            # ax_dna = plt.subplot(gs[1])
            # gff_file = file_path + chrom + ".gff"
            # # print(design) 

            # # Create the DNAplotlib renderer
            # design = dpl.load_design_from_gff(gff_file, chrom, region=cur_region)
            # dr = dpl.DNARenderer(scale=10.0)
            # part_renderers = dr.trace_part_renderers()
            # start, end = dr.renderDNA(ax_dna, design, part_renderers)
            # ax_dna.set_xlim(cur_region)
            # ax_dna.set_ylim([-8,8])
            # ax_dna.axis('off')

            # Alternative the DNAplotlib renderer
            # design = dpl.load_design_from_gff(gff_file, chrom, region=[cur_region[0]-5000,cur_region[1]+5000])
            # dr = dpl.DNARenderer(scale=10.0)
            # part_renderers = dr.SBOL_part_renderers()
            # start, end = dr.renderDNA(ax_dna, design, part_renderers)
            # ax_dna.set_xlim([start, end])
            # ax_dna.set_ylim([-15,15])
            # ax_dna.axis('off')

            # Update subplot spacing
            plt.subplots_adjust(hspace=0.2, left=0.15, right=0.85, top=0.85, bottom=0.15)
            
            # Save the figure
            fig.savefig(output_path + '\\' + chrom + '_' + str(cur_region[0]) + '.pdf', transparent=True)
            fig.savefig(output_path + '\\' + chrom + '_' + str(cur_region[0]) + '.png', dpi=300)

            # Clear the plotting cache
            plt.close('all')
            plt.clf()
            plt.cla()

def identify_promoter_TSS(self, RNA_seq_Profile, sample_name, Chromosome, Borders, alfa):
    
    profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
    RNAP_Flux_Profile = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile]
    min_value = min([x for x in RNAP_Flux_Profile if x > 7.2511e-10])
    diff_log = np.diff(np.log10([max(min_value, x) for x in RNAP_Flux_Profile]))
    diff_log_pos = [max(0, x) for x in diff_log]
    
    profile_rev = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
    RNAP_Flux_Profile_rev = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_rev]
    min_value_rev = min([x for x in RNAP_Flux_Profile_rev if x > 7.2511e-10])
    diff_log_rev = np.diff(np.log10([max(min_value_rev, x) for x in RNAP_Flux_Profile_rev]))
    diff_log_neg_rev = [min(0, x) for x in diff_log_rev]
    
    TSS_fwd = {sample_name:{Chromosome:[]}}
    TSS_rev = {sample_name:{Chromosome:[]}}
    
    for pos in range(len(diff_log_pos)):
        if diff_log_pos[pos] > 0.7:
            TSS_fwd[sample_name][Chromosome].append(Borders[0]+pos)
            
    for pos in range(len(diff_log_neg_rev)):
        if diff_log_neg_rev[pos] < -0.7:
            TSS_rev[sample_name][Chromosome].append(Borders[0]+pos)

    return diff_log_pos, diff_log_neg_rev, TSS_fwd, TSS_rev

def calculate_promoter_activity(self, RNA_seq_Profile, sample_name, Chromosome, Strand, TSS, X2, X1, alfa):

    Borders = [TSS-X2,
        TSS-X1,
        TSS+X1,
        TSS+X2]

    profile_upstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[0]:Borders[1]]
    UP = np.mean([alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_upstream])

    profile_downstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[2]:Borders[3]]
    DOWN = np.mean([alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_downstream])
            
    if Strand == 'fwd':
        promoter_activity = max(2e-7, DOWN - UP)
    elif Strand == 'rev':
        promoter_activity = max(2e-7, UP - DOWN)
        
    return promoter_activity

def identify_terminator_TTS(self, RNA_seq_Profile, sample_name, Chromosome, Borders, alfa, window_size):


    profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
    RNAP_Flux_Profile = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile]
    min_value = min([x for x in RNAP_Flux_Profile if x > 7.2511e-10])
    window_averaged_RNAP_Flux_Profile = [1 for kkk in range(window_size)] + [float(max(min_value,np.mean(RNAP_Flux_Profile[jjj-window_size:jjj]))) / max(min_value,np.mean(RNAP_Flux_Profile[jjj:jjj+window_size])) for jjj in range(window_size,len(RNAP_Flux_Profile)-window_size)] + [1 for kkk in range(window_size)]
    diff_log = np.log10(window_averaged_RNAP_Flux_Profile)
    diff_log_pos = [max(0, x) for x in diff_log]
    second_diff_log = np.diff([x if x > 0.7 else 0.0 for x in diff_log])
    
    TTS_fwd = {sample_name:{Chromosome:[]}}
    for wdw in range(1,len(second_diff_log)):
        if second_diff_log[wdw-1] > 0 and second_diff_log[wdw] <= 0:
            TTS_fwd[sample_name][Chromosome].append(Borders[0]+wdw)
    
    profile_rev = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
    RNAP_Flux_Profile_rev = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_rev]
    min_value_rev = min([x for x in RNAP_Flux_Profile_rev if x > 7.2511e-10])
    window_averaged_RNAP_Flux_Profile_rev = [1 for kkk in range(window_size)] + [float(max(min_value_rev,np.mean(RNAP_Flux_Profile_rev[jjj-window_size:jjj]))) / max(min_value_rev,np.mean(RNAP_Flux_Profile_rev[jjj:jjj+window_size])) for jjj in range(window_size,len(RNAP_Flux_Profile_rev)-window_size)] + [1 for kkk in range(window_size)]
    diff_log_rev = np.log10(window_averaged_RNAP_Flux_Profile_rev)
    diff_log_neg_rev = [min(0, x) for x in diff_log_rev]
    second_diff_log_rev = np.diff([x if x < -0.7 else 0.0 for x in diff_log_rev])
    
    TTS_rev = {sample_name:{Chromosome:[]}}
    for wdw in range(1,len(second_diff_log_rev)):
        if second_diff_log_rev[wdw-1] < 0 and second_diff_log_rev[wdw] >= 0:
            TTS_rev[sample_name][Chromosome].append(Borders[0]+wdw)
    
    TTS_window_fwd = {sample_name:{Chromosome:[]}}
    TTS_window_rev = {sample_name:{Chromosome:[]}}
    
    for pos in range(len(diff_log_pos)):
        if diff_log_pos[pos] > 0.7:
            TTS_window_fwd[sample_name][Chromosome].append(Borders[0]+pos)
            
    for pos in range(len(diff_log_neg_rev)):
        if diff_log_neg_rev[pos] < -0.7:
            TTS_window_rev[sample_name][Chromosome].append(Borders[0]+pos)

    return diff_log_pos, diff_log_neg_rev, TTS_window_fwd, TTS_window_rev, second_diff_log, second_diff_log_rev, TTS_fwd, TTS_rev

def calculate_terminator_strength(self, RNA_seq_Profile, sample_name, Chromosome, Strand, TTS, X2, X1, alfa):

		Borders = [TTS-X2,
			TTS-X1,
			TTS+X1,
			TTS+X2]

		profile_upstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[0]:Borders[1]]
		RNAP_Flux_Profile_up = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_upstream]
		UP = np.mean(RNAP_Flux_Profile_up)

		profile_downstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[2]:Borders[3]]
		RNAP_Flux_Profile_down = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_downstream]
		DOWN = np.mean(RNAP_Flux_Profile_down)
				
		if Strand == 'fwd':
			terminator_strength = float(UP)/DOWN
		elif Strand == 'rev':
			terminator_strength = float(DOWN)/UP
			
		return terminator_strength


# genomes = {'NZ_CP035288.1':2466502,'NZ_CP035289.1':20117,'NZ_CP035290.1':4439}
# # profile, stats = create_profile("Alignment_5A.SAM")
# # profile = load_profile("Alignment_5A"+"_profile.json",genomes)
# # norm_profile = normalize(profile,"Alignment_5A",'random.gff')
# profile = load_profile("Alignment_5A"+"_profile_norm.json",genomes)
# plot_regions_of_interest(profile,{'NZ_CP035288.1':[[1027400,1027500]]})


def read_gff(self, inputfiles, item):
    
    gene_dict = {aa:{} for aa in item}
    
    for inputfile in inputfiles:
        
        lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
        
        for line in lines:
            tokens = line.split('\t')
            genome, feature, start, end, strand = tokens[0], tokens[2], int(tokens[3]), int(tokens[4]), tokens[6]
            name = (tokens[-1].split('Name='))[-1]
            gene_dict[genome][name] = {'start': start, 'end': end, 'strand': strand, 'feature': feature}
            
    return gene_dict

def map_raw_FASTQ_files(sample_name,RAW_FASTQ_filenames,Fasta_filename,ebwt_base):

        Trimmed_RAW_FASTQ_filename = default_prefix + 'Trimmed_' + sample_name
        Unaligned_FASTQ_filename = default_prefix + 'Unaligned_' + sample_name
        
        cmd_index = 'AdapterRemoval --file1 ' + RAW_FASTQ_filenames[0] + ' --file2 ' + RAW_FASTQ_filenames[1] + ' --basename ' + sample_name
        # 'cutadapt -a CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG —-discard-untrimmed -m 10 -o ' +  + ' ' + sample_name
        print("Removing linker: "+cmd_index)
        status0 = subprocess.call(cmd_index, shell=True)
        
        #Linker removal
        #cutadapt (download site: https://cutadapt.readthedocs.org/en/stable/)
        cmd_index = 'cutadapt -a CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG —-discard-untrimmed -m 10 -o ' + Trimmed_RAW_FASTQ_filename + ' ' + RAW_FASTQ_filename[sample_name]
        print("Removing linker: "+cmd_index)
        status0 = subprocess.call(cmd_index, shell=True)

        # ----------------------------------- #
        
        ## Make the indexes
        #cmd_index = 'bowtie-build' + ' -p ' + self.NonCoding_Genes_Fasta_filename[sample_name] + ' ' + self.NonCoding_Genes_INDEXED_Fasta_Filename[sample_name]
        #print("Making index Non-coding genes: "+cmd_index)
        #status1 = subprocess.call(cmd_index, shell=True)
        
        # # Make the indexes
        # cmd_index = 'C://Users//Hamid//bowtie//bowtie2-build ' + default_prefix+Fasta_filename + ' ' + default_prefix+ebwt_base
        # print("Making index coding genes: "+cmd_index)
        # status2 = subprocess.call(cmd_index, shell=True)
        
        # ----------------------------------- #
        
        ## Perform the mapping
        #cmd_mapping = 'bowtie' + ' —-best -t —-un ' + self.Coding_Gene_FASTQ_filename[sample_name] + ' ' + self.NonCoding_Genes_INDEXED_Fasta_Filename[sample_name] + ' ' + self.Trimmed_RAW_FASTQ_filename[sample_name] + ' ' + self.NonCoding_Genes_MAP_output_filename[sample_name]
        #print("Mapping Reads Bowtie: "+cmd_mapping)
        #status3 = subprocess.call(cmd_mapping, shell=True)
        
        # Perform the mapping
        cmd_mapping = 'bowtie2 --local —S' + default_prefix+sample_name + '.sam --best -x' + default_prefix+'Staphepi' + ' -1 ' + default_prefix+RAW_FASTQ_filenames[0] + ' -2 ' + default_prefix+RAW_FASTQ_filenames[1]
        print("Mapping Reads Bowtie: "+cmd_mapping)
        status4 = subprocess.call(cmd_mapping, shell=True)
        
        # Perform the mapping
        # cmd_mapping = 'bowtie' + ' —t -v1 -m2 -k1 —-un ' + Unaligned_FASTQ_filename + ' ' + self.Coding_Genes_INDEXED_Fasta_Filename[sample_name] + ' ' + Trimmed_RAW_FASTQ_filename + ' ' + Coding_Genes_MAP_output_filename
        # print("Mapping Reads Bowtie: "+cmd_mapping)
        # status4 = subprocess.call(cmd_mapping, shell=True)

map_raw_FASTQ_files('15B',['5B_R_1.fastq','5B_R_2.fastq'],'GCF_006094375.1_ASM609437v1_genomic.fa','GCF_006094375')
