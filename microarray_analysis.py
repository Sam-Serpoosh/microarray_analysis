from read_csv import create_genes_values_dictionary
from math import sqrt
from math import fabs
from random import shuffle
import numpy as np
import operator

NUM_OF_DISEASE = 3
NUM_OF_PERMUTATIONS = 10

def find_discriminant_genes(patients_info_file):
  genes_values = create_genes_values_dictionary(patients_info_file)
  genes_t_stats = { k: [] for k in genes_values.keys() }
  genes_d_scores = {}
  for gene_name, gene_values in genes_values.items():
    genes_t_stats[gene_name].append(calculate_t_stat(gene_values))
    randomly_distribute_and_calculate_t_stat(genes_t_stats, 
        gene_name, gene_values)

  for gene_name, t_stats in genes_t_stats.items():
    genes_d_scores[gene_name] = calculate_d_score(t_stats)

  sorted_d_scores = sorted(genes_d_scores.items(), 
      key=operator.itemgetter(1), reverse=True)
  print(sorted_d_scores)

def randomly_distribute_and_calculate_t_stat(genes_t_stats, 
    gene_name, gene_values):
  for i in range(0, NUM_OF_PERMUTATIONS):
    shuffle(gene_values)
    t_stat = calculate_t_stat(gene_values)
    genes_t_stats[gene_name].append(t_stat)

def calculate_t_stat(values):
  disease_part = values[0:NUM_OF_DISEASE]
  normal_part = values[NUM_OF_DISEASE:len(values)]
  disease_values = [num for num in disease_part if num != -99]
  normal_values = [num for num in normal_part if num != -99]

  return t_stat(disease_values, normal_values)


def t_stat(values1, values2):
  np_values1 = np.array(values1)
  np_values2 = np.array(values2)
  mean1 = np_values1.mean()
  mean2 = np_values2.mean()
  std1 = np.std(np_values1, ddof=1, dtype=np.float32)
  std2 = np.std(np_values2, ddof=1, dtype=np.float32)

  numerator = mean1 - mean2
  denominator = sqrt((std1 ** 2 / len(np_values1)) + 
      (std2 ** 2 / len(np_values2)))

  return numerator / denominator

def calculate_d_score(t_stats):
  distribution_t_stats = np.array(t_stats[1:len(t_stats)])
  distribution_mean = distribution_t_stats.mean()
  distribution_std = np.std(distribution_t_stats, ddof=1, dtype=np.float32)

  return fabs(t_stats[0] - distribution_mean) / distribution_std

find_discriminant_genes("working_sample.csv")
