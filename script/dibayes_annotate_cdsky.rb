#!/usr/bin/env ruby
# 
# dibayes_annotate_cdsky.rb
# A script that read dibayes result and annotate each SNP
# based on KY decision file

require 'optparse'
require 'bio'
require 'yaml'

options = {}

optparse = OptionParser.new do |opts|

  opts.banner = "Usage: #{File.basename($0)} [options] dibayes.tab ky_decision_file"

  opts.on('-c', '--conf FILE', 'config file - required') do |v|
    options[:conf] = v
  end

  opts.on('-o', '--out FILE', "output file - optional") do |v|
    options[:out] = v
  end

  opts.on('-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end

end

optparse.parse!

if ARGV.size == 2
  $dibayes_tab = ARGV[0]
  $ky_annot_file = ARGV[1]
else
  raise "\nArgument Error: 'dibayes.tab' and 'ky_decision_file' should be specified.\n#{optparse}"
end

## load conf file
unless options[:conf]
  raise "\nArgument Error: conf file required.\n#{optparse}\n"
end
conf = YAML.load(File.open(options[:conf]).read)
  




$outf = (conf[:out] || "#{File.basename($dibayes_tab)}.ann.txt")
o = File.open($outf, "w")

STDERR.puts "Start processing ..."
STDERR.puts "Annotated file is being saved as '#{$outf}'"

## Added columns are
# 0: gene model ID
# 1: strand
# 2: codon change pattern
# 3: aa change pattern
# 4: position in aa
# 5: position in na
# 6: change type
# 7: TAIR9 description

NEW_COLNAMES = %w{
gene_model
strand
codon_change
aa_change
position_aa
position_nuc
change_type
description
}

def seqid2chr(seqid)
  seqid.sub(/^Chr/, '')
end

## Read Yamaguchi's annotation table

# 0:    chhr-position	1_100997
# 1:	chr-position-seq	1_100997_C
# 2	position	100997
# 3	chr	1
# 4	coding_strand	+
# 5	locus_protein_code	AT1G01240.1
# 6	base_seq	C
# 7	ref_seq	G
# 8	sample_seq_num	315
# 9	coding_nuc_change	TAC=>TAX
# 10	ref_AA	Y
# 11	coding_protein_num	105
# 12	AA_seq(X=>ACGT)	*Y*Y
# 13	coding_protein_length	331
# 14	ref_AA	Y
# 15	sample_coding_AA	*
# 16	dicision	mutable

data = Hash.new
## key = [chr, position]
## value = array of all elements 0-16

File.open($ky_annot_file).each do |l|
  a = l.chomp.split(/\t/)
  chr = a[3]
  pos = a[2].to_i
  key = [chr, pos]
  unless data.has_key?(key)
    data[key] = []
  end
  data[[chr, pos, ]] << a
end

## Load TAIR9 description
tair9desc = Hash.new
File.open(conf["tair9_desc_file"]).each do |l|
  a = l.chomp.split(/\t/)
  name = a[0]
  tair9desc[name] = {
    'short_desc' => a[2],
    'summary' => a[3]
  }
                                         
end

## Process DiBayes table

$> = o

## print header
puts "#=== SNP Annotation by santenna ==="
puts "#"
puts "# format:  santenna_cds 1.0"
puts "# module:  #{conf['module']}"
puts "# source1: #{$dibayes_tab}"
puts "# source2: #{$ky_annot_file}"
puts "# genome:  #{conf['ref_genome_name']}"
puts "# script:  #{File.basename($0)}"
puts "# operator:#{ENV['USER']}"
puts "# date:    #{Time.now}"
puts "#"
puts "# copyright: NIBB Core Research Facilities"
puts "#"

File.open($dibayes_tab).each_with_index do |l, i|
  a = l.chomp.split(/\t/)

  if i == 0
    col_chosen = %w{chromosome  position  reference  genotype  score
                    coverage 
                    refAlleleCounts novelAlleleCounts refAlleleMeanQV novelAlleleMeanQV
                    homo_hetero }
    puts "#" + col_chosen.concat(NEW_COLNAMES).join("\t")
    next
  end

  seqid = a[0]
  pos = a[3].to_i
  chr = seqid2chr(seqid)
  kydata = data[[chr, pos]]
  if kydata
    kydata.each do |d|
      new_allele = d[7].delete(d[6])
      strand = d[4]

      codon_ref, codon_new = d[9].downcase.split(/\->/)
      new_allele2 = Bio::Sequence::NA.new(new_allele)
      if strand == "-"
        new_allele2 = new_allele2.complement
      end

      codon_new = codon_new.sub(/x/, new_allele2.to_s.upcase)

      aa_ref = Bio::Sequence::NA.new(codon_ref).translate
      aa_new = Bio::Sequence::NA.new(codon_new).translate

      
#      [codon_ref, codon_new, aa_ref, aa_new]

      aa_len = d[13].to_i
      aa_pos = d[11].to_i
      

      nuc_len = aa_len * 3  + 3 # stop codon included
#      nuc_pos = (aa_pos - 1) * 3 + codon_new.index(/[ATGC]/) + 1
      nuc_pos = d[8].to_i

#p      [aa_len, aa_pos, nuc_len, nuc_pos]

      change_type = (d[16] == "mutable" ? "NonSynonymous" : "Synonymous")

      gene = d[5]
      genotype = d[7]
      if genotype.size == 1
        genotype = genotype * 2
      end

       out = []
       out << [a[0], #seqid
               a[3].to_i, #start == position
               d[6], #reference
               genotype,
               a[5], #score
               a[10].to_i, #coverage
               a[11].to_i, #refAlleleCounts
               a[14].to_i, #novelAlleleCounts
               a[13].to_i, #refAlleleMeanQV
               a[16].to_i, #novelAlleleMeanQV
               (a[19] == "1" ? "het" : "homo") #hetero or homo
              ]

       out << [gene, #gene model ID
               d[4], #strand
               "#{codon_ref}=>#{codon_new}" , #codon change pattern
               "#{aa_ref}=>#{aa_new}", #aa change pattern
               "#{aa_pos}/#{aa_len}", #position in aa
               "#{nuc_pos}/#{nuc_len}", #position in nuc
               change_type
              ]

       out << [tair9desc[gene]['short_desc']] 

       puts out.flatten.join("\t")
     end
   end

end

STDERR.puts "done"
