import requests
import csv
import StringIO
import urllib
import os
import sys
import pprint

class CaseSet(object):
	genetic_profile_ids = ['gbm_tcga_gistic', 'gbm_tcga_mutations']
	cases = {}
	gene_meta_data = {}
	meta_data = {}
	summary = {}

	def __init__(self, gene_list):
		"""Initialized the case-set object by adding data for the GBM genetic profiles data for for each gene in the gene_list parameter"""
		for genetic_profile_id in self.genetic_profile_ids:
			self.add_genetic_profile(genetic_profile_id, gene_list, useLocalFile=False)

	def add_genetic_profile(self, genetic_profile_id, gene_list, useLocalFile=False):
		"""Adds a genetic profile for every gene provided."""
		gene_alterations = {'meta_data': [], 'gene_meta_data': {}, 'cases': {} }
		gene_meta_data_fields = ['COMMON', 'GENE_ID']
		
		if type(gene_list) == list():
			gene_list = ','.join(gene_list)
		
		if useLocalFile:
			file_obj = open(self.get_profile_data_file_name(genetic_profile_id, gene_list), 'r')
		else:
			url = self.get_profile_data_URL(genetic_profile_id, gene_list)
			# create a object with same interface as open file to allow for the use of file-based parsing methods
			file_obj = self.get_profile_data_API_response(url)

		# create an array to store the position of each line the file, to allow for multiple rows of meta-data to proceed the column headers
		header_row_begin_position = 0
		for line in iter(file_obj.readline, ''):
			# check if the line begins with "#" as this is the convention for providing meta-data in the text returned by the API
			if line[0] == '#':
				header_row_begin_position = file_obj.tell()
			else:
				break
		# get the position before the header row (the file curesor will be at the end of the end of the header row when the loop exits)
		file_obj.seek(header_row_begin_position)
		reader = csv.DictReader(file_obj, delimiter='\t')
		for row in reader:
			# get the common abbreviation for the gene represented in the row of data
			gene = row['COMMON']
			for field in row:
				# populate the meta data for gene represented in this row
				if field in gene_meta_data_fields:
					if not (gene in self.gene_meta_data):
						self.gene_meta_data[gene] = {}
					self.gene_meta_data[gene][field] = row[field]
				else:
					if not (field in self.cases):
						self.cases[field] = {}
					if not (gene in self.cases[field]):
						self.cases[field][gene] = {}
					self.cases[field][gene][genetic_profile_id] = row[field]
		file_obj.close()

	def display_case_set_summary(self):
		""""Displays a summary of all the analysis done on the GBM genetic profile data available for all specfied genes."""
		display = [{'analysis': 'profile_gbm_tcga_mutations_summary', 'text': '{gene} is mutated in {result_text}% of all cases.', 'summary_function': self.get_mutation_percent},
					{'analysis': 'profile_gbm_tcga_gistic_summary',  'text': '{gene} is copy number altered in {result_text}% of all cases.', 'summary_function': self.get_copy_number_altered_percent},
					{'analysis': 'aggregate', 'text': 'Total % of cases where {gene} is altered by either mutation or copy number alteration: {result_text}% of all cases.', 'summary_function': self.get_all_alterated_percent}]
		
		self.set_summary()
		tokens = {}
		for gene in self.summary:
			tokens = {'gene': gene}
			for metric in display:
				result_percent= metric['summary_function'](gene)
				tokens['result_text'] = '{0:.{1}f}'.format(result_percent*100, 0)
				print metric['text'].format(**tokens)

	def get_mutation_percent(self, gene):
		"""Returns the perceantage (as  decimal) of all cases where mutations were detected for a specific gene (includes cases with multiple alterations)."""
		return float(self.summary[gene]['mutated_case_count']) / float(self.summary[gene]['total_case_count'])

	def get_copy_number_altered_percent(self, gene):
		"""Returns the perceantage (as  decimal) of all cases where copy number alterations were detected for a specific gene (includes cases with multiple alterations)."""		
		return float(self.summary[gene]['copy_number_alterated_case_count']) / float(self.summary[gene]['total_case_count'])

	def get_all_alterated_percent(self, gene):
		"""Returns the perceantage (as  decimal) of all cases where alterations (mutation or copy number alterations) where detected for a specific gene (includes cases with multiple alterations)."""
		return float(self.summary[gene]['copy_number_alterated_case_count'] + self.summary[gene]['mutated_case_count'] - self.summary[gene]['multiple_alterations_case_count']) / float(self.summary[gene]['total_case_count'])

	def is_mutated(self, case, gene):
		"""Returns a boolean indicating if the a specific gene is mutated for the case in question."""
		genetic_profile_id = 'gbm_tcga_mutations'
		profile_result = case[gene][genetic_profile_id] 
		if profile_result == 'NaN' or str(profile_result) == '0':
			return False
		else:
			return True

	def is_copy_number_alterated(self, case, gene):
		"""Returns a boolean indicating if the a specific gene has copy number alteration for the case in question."""		
		genetic_profile_id = 'gbm_tcga_gistic'
		profile_result = case[gene][genetic_profile_id]
		# @TODO: determine best way to handle 'NA' values
		if int(profile_result) not in [0, 1, -1]:
			return True
		else:
			return False

	def set_summary(self):
		"""Summarizes the alterations (mutations or copy number alterations) observerd for all cases in the case set,
		 looking at all genes specified at the creation of CaseSet objet."""
		mutated_case_count = 0
		copy_number_alterated_case_count = 0
		multiple_alterations_case_count = 0		
		self.total_case_count = len(self.cases)
		for gene in self.gene_meta_data:
			for case in self.cases:
				if self.is_copy_number_alterated(self.cases[case], gene):
					copy_number_alterated_case_count += 1
				if self.is_mutated(self.cases[case], gene):
					mutated_case_count += 1
				if self.is_copy_number_alterated(self.cases[case], gene) and self.is_mutated(self.cases[case], gene):
					multiple_alterations_case_count += 1
			
			mutated_case_percent = float(mutated_case_count) / float(self.total_case_count)
			copy_number_alterated_percent = float(copy_number_alterated_case_count) / float(self.total_case_count)
			all_alterated_percent = float(copy_number_alterated_case_count + mutated_case_count - multiple_alterations_case_count) / float(self.total_case_count)
			self.summary[gene] = {'total_case_count': self.total_case_count,
										'mutated_case_count': mutated_case_count,
										'copy_number_alterated_case_count': copy_number_alterated_case_count,
										'multiple_alterations_case_count': multiple_alterations_case_count,
										'all_alterated_percent': all_alterated_percent}

	def get_request_payload(self, genetic_profile_id, gene_list):
		"""Returns the correctly formatted paylod (query params) for making request the cBioPortal web API."""
		payload = {}
		payload['cmd'] = 'getProfileData'
		payload['id_type'] = 'gene_symbol'
		payload['case_set_id'] = 'gbm_tcga_cnaseq'
		payload['genetic_profile_id'] = genetic_profile_id
		payload['gene_list'] = gene_list
		return payload

	def get_profile_data_API_response(self, url):
		"""Returns an object for interacting with data povided by the API and exposing the same interface as an open file object."""
		response = requests.get(url)
		file_obj = StringIO.StringIO(response.text)
		return file_obj

	def get_profile_data_file_name(self, genetic_profile_id, gene_list):
		"""This function returns the file_path to data cahced from the web API.
		Cached files simply use the URL-encoded querystring of a call to the web API for the file name (with .txt as an extension)"""
		file_dir = '../data_for_cbioportal_cmd_util/'
		payload = self.get_request_payload(genetic_profile_id, gene_list)
		query = 'cmd={cmd}&genetic_profile_id={genetic_profile_id}&id_type={id_type}&gene_list={gene_list}&case_set_id={case_set_id}'.format(**payload)
		file_name = file_dir + urllib.quote(query) + '.txt'
		return file_name

	def get_profile_data_URL(self, genetic_profile_id, gene_list):
		"""Returns a the URL to be used to query the web API for the relevant data."""
		payload = self.get_request_payload(genetic_profile_id, gene_list)	
		return 'http://www.cbioportal.org/webservice.do?cmd={cmd}&genetic_profile_id={genetic_profile_id}&id_type={id_type}&gene_list={gene_list}&case_set_id={case_set_id}'.format(**payload)

def main():
	"""Executes whenever the script is executed directly (not imported). Handles arguments passed from the command line."""
	if ( len(sys.argv) == 1 ) or ( len(sys.argv) > 4 ):
		print('usage: gbm_summarize.py gene1 [gene2, gene3]')
	else:
		# get all arguments except the argv[0] (the script file), add them to a comma-delimited string
		gene_list = ','.join(sys.argv[1:])
		case_set = CaseSet(gene_list)
		case_set.display_case_set_summary()

if __name__ == '__main__':
	main()