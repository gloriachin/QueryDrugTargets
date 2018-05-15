#!/usr/bin/env python
import xml.etree.ElementTree as ET
import sys
import codecs
#reload(sys)

fin = open('/ref/UniprotId_GeneName_RefseqId.TXT','r')
dic_AC_GN={}
dic_GN_AC = {}
for line in fin.readlines():
    words = line.split("\t")
    leng = len(words)
    if leng > 2:
      AC = words[0].strip(';')
      GN_pre = words[1].split(';')[0]
      if '{' in GN_pre:
        GN = GN_pre.split('{')[0].strip(' ')
      else:
        GN = GN_pre
      dic_AC_GN[AC] =GN
      dic_GN_AC[GN] = AC
fin.close()

fin = open("/Data/Protein/expression.txt",'r')
dic_ac_pro_name={}
dic_GN_pro_name={}

for line in fin.readlines():
  words = line.strip().split('\t')
  ac = words[1]
  pro_name = words[4]
  dic_ac_pro_name[ac]=pro_name
  if ac in dic_AC_GN:
    gn = dic_AC_GN[ac]
    pro_name = words[4]
    dic_GN_pro_name[gn] = pro_name
fin.close()



dic_TTD_approved={}
fin = open("/ref/TTD/TTD_uniprot_success_new.txt")
for line in fin.readlines():
  key = line.strip().split('\t')[-1]
  dic_TTD_approved[key]=''
fin.close()

dic_TTD_clinical={}
fin = open("/Users/gloria/Documents/Project/coding_projects/TTD/TTD_uniprot_clinical.txt")
for line in fin.readlines():
  key = line.strip().split('\t')[-1]
  dic_TTD_clinical[key]=''
fin.close()

dic_TTD_research={}
fin = open("ref/TTD/TTD_uniprot_research.txt")
for line in fin.readlines():
  key = line.strip().split('\t')[-1]
  dic_TTD_research[key]=''
fin.close()

dic_query={}
fin = open(sys.argv[1],'r')
for line in fin.readlines():
  query = line.strip()
  if query in dic_AC_GN:
    dic_query[query]=''
  elif query in dic_GN_AC:
    query = dic_GN_AC[query]
    dic_query[query] = ''
fin.close()



print("start parse....")
#sys.setdefaultendcoding('utf-8')
#inF = '/Users/gloria/Documents/Project/coding_projects/Drugbank/Drugbank_full_database.xml'
inF = '/ref/Drugbank/full_database_2017.11.3.xml'
root = ET.parse(inF)
namespaces = {'ns': 'http://www.drugbank.ca'}
print("root.findall....")

druglist = root.findall("./ns:drug",namespaces)
id = 0
fou = codecs.open(sys.argv[1]+'.drug_all.csv','w', encoding='utf-8')
for drug in druglist:
  id = id + 1

  name = drug.find("./ns:name", namespaces)
  bid = drug.find("./ns:drugbank-id[@primary]", namespaces)
  cas = drug.find("./ns:cas-number", namespaces)
  targets = drug.findall("./ns:targets/ns:target", namespaces)
  group = drug.find("./ns:groups/ns:group", namespaces)
  indication = drug.find("./ns:indication", namespaces)
  pharmacodynamics = drug.find("./ns:pharmacodynamics", namespaces)
  toxicity = drug.find("./ns:toxicity", namespaces)
  externalidentifier = drug.find("./ns:external-identifiers/ns:external-identifier", namespaces)

  gene_name_list = drug.findall("./ns:targets/ns:target/ns:polypeptide[@id]",namespaces)
  tlist = []

  for target in targets:
    t = {}
    t['name'] = target.find("./ns:name", namespaces).text
    alist = []
    for action in target.findall("./ns:actions/ns:action", namespaces):
      alist.append(action.text)
    t['actions'] = alist
    tlist.append(t)

  #if group.text == "approved":
  if gene_name_list != None:
    for gene_name in gene_name_list:
      if gene_name.attrib['id'] in dic_query:
        fou.write(name.text+'\t'+action.text+'\t'+gene_name.attrib['id']+'\t'+dic_AC_GN[gene_name.attrib['id']]+'\t'+group.text+'\t'+dic_ac_pro_name[gene_name.attrib['id']])
        if gene_name.attrib['id'] in dic_TTD_approved:
          fou.write('\t'+'TTD_approved\n')
        elif gene_name.attrib['id'] in dic_TTD_clinical:
          fou.write('\t'+'TTD_clinical\n')
        elif gene_name.attrib['id'] in dic_TTD_research:
          fou.write('\t' + 'TTD_research\n')
        else:
          fou.write('\n')
fou.close()
