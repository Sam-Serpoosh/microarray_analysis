from read_csv import create_genes_values_dictionary
from math import sqrt
from math import fabs
from math import isnan
from random import shuffle
from mpi4py import MPI
import numpy as np
import operator

NOT_EXIST = -99
NUM_OF_DISEASE = 8
NUM_OF_PERMUTATIONS = 1000

def find_discriminant_genes(patients_info_file):
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  genes_values_chunks = chunk_genes_values_for_scattering(
          patients_info_file, comm, rank, size)
  genes_d_scores = calculate_all_genes_d_scores(comm, genes_values_chunks)
  if rank == 0:
      sorted_d_scores = sort_genes_d_scores(genes_d_scores)
      write_scores_to_file("discriminant_genes.txt", sorted_d_scores)


def chunk_genes_values_for_scattering(patients_info_file, comm, rank, size):
  if rank == 0: # Master
    genes_values = create_genes_values_dictionary(patients_info_file)
    genes_values_chunks = [{} for _ in range(size)]
    index = 0
    for gene_name, gene_values in genes_values.items():
      genes_values_chunks[index % size][gene_name] = gene_values
      index += 1
    return genes_values_chunks
  else:
    return None

def calculate_all_genes_d_scores(comm, genes_values_chunks):
  genes_chunk = comm.scatter(genes_values_chunks, root=0)
  genes_t_stats = dict((k, []) for k in genes_chunk.keys())
  for gene_name, gene_values in genes_chunk.items():
    genes_t_stats[gene_name].append(calculate_t_stat(gene_values))
    randomly_distribute_and_calculate_t_stat(genes_t_stats,
        gene_name, gene_values)

  genes_d_scores = {}
  for gene_name, t_stats in genes_t_stats.items():
    genes_d_scores[gene_name] = calculate_d_score(t_stats)

  return comm.gather(genes_d_scores, root=0)

def randomly_distribute_and_calculate_t_stat(genes_t_stats,
    gene_name, gene_values):
  for i in range(0, NUM_OF_PERMUTATIONS):
    shuffle(gene_values)
    t_stat = calculate_t_stat(gene_values)
    genes_t_stats[gene_name].append(t_stat)

def calculate_t_stat(values):
  disease_part = values[0:NUM_OF_DISEASE]
  normal_part = values[NUM_OF_DISEASE:len(values)]
  disease_values = [num for num in disease_part if num != NOT_EXIST]
  normal_values = [num for num in normal_part if num != NOT_EXIST]

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

  score = fabs(t_stats[0] - distribution_mean) / distribution_std
  if isnan(score):
    return NOT_EXIST
  return score

def sort_genes_d_scores(genes_d_scores):
    all_d_scores = {}
    for d_scores in genes_d_scores:
      for gene_name, score in d_scores.items():
        all_d_scores[gene_name] = score

    return sorted(all_d_scores.items(), key=lambda x:x[1], reverse=True)

def write_scores_to_file(output_file, scores):
    f = open(output_file, "w")
    for score in scores:
      f.write(str(score) + "\n")
    f.close()

find_discriminant_genes("NCI-60.csv")
