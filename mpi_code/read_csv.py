import csv
import io

def create_genes_values_dictionary(filename):
  genes_values = {}
  with io.open(filename, "r", newline="") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=",")
    counter = 0
    for row in csv_reader:
      if counter == 0:
        counter += 1
        continue
      add_gene_name_and_values(genes_values, row)

  return genes_values

def add_gene_name_and_values(genes_values, gene_info):
  gene_name = gene_info.pop(0).strip()
  gene_info.pop(0) # get rid of the Unique Identifier
  for index, item in enumerate(gene_info):
    if (item.strip() == ""):
      gene_info[index] = "-99"

  values = [float(num.strip()) for num in gene_info]
  genes_values[gene_name] = values
