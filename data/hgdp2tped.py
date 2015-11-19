# Python script to convert HGDP genotype data to the transposed PLINK
# (.tped) format. Modified slightly from the script originally written
# by Prem Gopalan, Wei Hao, David Blei and John Storey as part of the
# TeraStructure repository on Github:
#
#   http://github.com/premgopalan/terastructure
#

# Read the SampleInformation.txt file. Note that the last line of the
# file SampleInformation.txt (giving column totals) should have been
# removed.
f = open("SampleInformation.txt","r")

# Get the list of ids in the "H952" set.
SampleInformationHeader = f.readline().strip().split()
i1 = SampleInformationHeader.index("IncludedInDataset952[No1stOr2ndDegreeRelatives]?")
i2 = SampleInformationHeader.index("PopulationName")
indiv = {}

for x in f:
  y = x.strip().split()
  if y[i1] == '1':
    tmp = '000000000000' + y[0]
    tmp = tmp[-5:]
    indiv['HGDP' + tmp] = y[i2]

f.close()

# These are the ids in H952 but not in website data.
del indiv['HGDP00987']
del indiv['HGDP00453']
del indiv['HGDP00452']
del indiv['HGDP01219']
del indiv['HGDP00247']
del indiv['HGDP00248']
del indiv['HGDP01149']
del indiv['HGDP00660']
del indiv['HGDP01344']
del indiv['HGDP01233']
del indiv['HGDP00754']
del indiv['HGDP00944']

geno_in  = open("HGDP_FinalReport_Forward.txt","r")
map_in   = open("HGDP_Map.txt","r")
geno_out = open("hgdp.tped","w")

hgdp_indiv = geno_in.readline()
hgdp_indiv = hgdp_indiv.strip()
hgdp_indiv = hgdp_indiv.split()
valid      = [hgdp_indiv.index(x) for x in indiv.keys()]

# Convert the genotype data to the transposed PLINK format.
for x in geno_in:
  snpinfo = map_in.readline().split()
  if snpinfo[1] not in map(str,range(1,23)):
    continue
  x = x.strip()
  x = x.split()
  x = x[1:]
  geno = [x[y] for y in valid]
  out_tped = [snpinfo[1],snpinfo[0],"1",snpinfo[2]]
  missing_count = 0.0
  allele1_count = 0.0
  tot_allele    = 0.0
  homozygous1   = ''
  for y in xrange(0,len(geno)):
    if geno[y] == '--':
      out_tped.append("0 0")
      missing_count = missing_count + 1
    else:
      out_tped.append(geno[y][0] + ' ' + geno[y][1])
      tot_allele = tot_allele + 2.0
      if geno[y][1] != geno[y][0]:
        allele1_count = allele1_count + 1.0
      elif homozygous1 == '':
        homozygous1   = geno[y]
        allele1_count = allele1_count + 2.0
      elif homozygous1 == geno[y]:
        allele1_count = allele1_count + 2.0

  # At most 5% missing genotypes per SNP allowed.
  if missing_count <= 47: 
    if allele1_count/tot_allele >= 0.01 or allele1_count/tot_allele <= 0.99:
      geno_out.write(" ".join(out_tped) + "\n")
    else:
      print("error---this shouldn't happen in HGDP data")

geno_in.close()
map_in.close()
geno_out.close()

# Write out the .tfam file.
f = open("hgdp.tfam","w")

for x in valid:
  out_tfam = [str(x),hgdp_indiv[x],'0','0','1','0']
  f.write(" ".join(out_tfam) + "\n")
f.close()
