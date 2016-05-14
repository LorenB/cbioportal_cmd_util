import requests
import csv
import StringIO
import urllib
import os
import sys
import pprint

genetic_profile_lookup = {
	'mutations': {'genetic_profile_id': 'gbm_tcga_mutations', 'display_text': '{COMMON} is mutated in {percent_display} of all cases.'},
	'copy_number_alterations': {'genetic_profile_id': 'gbm_tcga_gistic',  'display_text': '{COMMON} is copy number altered in {percent_display} of all cases.'},
	'aggregate': {'display_text': 'Total % of cases where {COMMON} is altered by either mutation or copy number alteration: {percent_display} of all cases.'}
}

def display_gene_alterations_summary_descriptions(gene_list):
	metrics = ['mutations', 'copy_number_alterations', 'aggregate']
	for gene in gene_list:
		gene_alterations_summary_descriptions = get_gene_alterations_summary_descriptions(gene)
		for metric in metrics:
			print gene_alterations_summary_descriptions[metric]		

def get_gene_alterations_summary_descriptions(gene):
	gene_alterations_data = get_all_gene_alterations_data(gene)
	gene_alterations_summary = get_all_gene_alterations_summary(gene_alterations_data)
	gene_alterations_summary_descriptions = {}
	display_data = {}
	for alteration_type in gene_alterations_summary['data']:
		display_data = gene_alterations_summary['gene_meta_data'].copy() 
		display_data['percent_display'] = gene_alterations_summary['data'][alteration_type]['percent_display']
		gene_alterations_summary_descriptions[alteration_type] = genetic_profile_lookup[alteration_type]['display_text'].format(**display_data) 
	return gene_alterations_summary_descriptions

def get_all_gene_alterations_data(gene):
	gene_alterations_data = {} 
	gene_alterations_data['mutations'] = get_gene_alterations(genetic_profile_lookup['mutations']['genetic_profile_id'], gene)
	gene_alterations_data['copy_number_alterations'] = get_gene_alterations(genetic_profile_lookup['copy_number_alterations']['genetic_profile_id'], gene)
	return gene_alterations_data

def get_all_gene_alterations_summary(gene_alterations_data):
	gene_alterations_summary = {'data': {}, 'gene_meta_data': {} }
	gene_alterations_summary['gene_meta_data'] = gene_alterations_data['mutations']['gene_meta_data']
	gene_alterations_summary['data']['mutations'] = get_gene_mutation_data_summary(gene_alterations_data['mutations']['cases'])
	gene_alterations_summary['data']['copy_number_alterations'] = get_gene_copy_number_alteration_data_summary(gene_alterations_data['copy_number_alterations']['cases'])	
	gene_alterations_summary['data']['aggregate'] = {}
	all_altered_cases = set(gene_alterations_summary['data']['mutations']['altered_case_ids']) | set(gene_alterations_summary['data']['copy_number_alterations']['altered_case_ids'])
	all_altered_cases_count = len(all_altered_cases)
	total_cases = len(gene_alterations_data['mutations']['cases'])
	gene_alterations_summary['data']['aggregate']['percent_decimal'] = float(all_altered_cases_count)/float(total_cases)
	gene_alterations_summary['data']['aggregate']['percent_display'] = str( int( round(gene_alterations_summary['data']['aggregate']['percent_decimal']*100, 0) ) ) + r'%'
	
	return gene_alterations_summary

def get_gene_alterations(genetic_profile_id, gene, useLocalFile=True):
	'''Returns one genetic profile for ever gene provided.'''
	gene_alterations = {'meta_data': [], 'gene_meta_data': {}, 'cases': {} }
	gene_meta_data_fields = ['COMMON', 'GENE_ID']
	url = get_profile_data_URL(genetic_profile_id, gene)
	
	if useLocalFile:
		file_obj = get_profile_data_from_file(get_profile_data_file_name(genetic_profile_id, gene))
	else:
		# create a object with same interface as open file to allow for the use of file-based parsing methods
		file_obj = get_profile_data_API_response(url)

	# create an array to store the position of each line the file, to allow for multiple rows of meta-data to proceed the column headers
	header_row_begin_position = 0
	for line in iter(file_obj.readline, ''):
		# check if the line begins with "#" as this is the convention for providing meta-data in the text returned by the API
		if line[0] == '#':
			gene_alterations['meta_data'].append(line)
			header_row_begin_position = file_obj.tell()
		else:
			break
	# get the position before the header row (the file curesor will be at the end of the end of the header row when the loop exits)
	file_obj.seek(header_row_begin_position)
	reader = csv.DictReader(file_obj, delimiter='\t')
	row = reader.next()

	for field in row:
		if field in gene_meta_data_fields:
			gene_alterations['gene_meta_data'][field] = row[field]
		else:
			gene_alterations['cases'][field] = row[field]
	file_obj.close()
	return gene_alterations

def get_request_payload(genetic_profile_id, gene_list):
	payload = {}
	payload['cmd'] = 'getProfileData'
	payload['id_type'] = 'gene_symbol'
	payload['case_set_id'] = 'gbm_tcga_cnaseq'
	payload['genetic_profile_id'] = genetic_profile_id
	payload['gene_list'] = gene_list
	return payload

def get_gene_copy_number_alteration_data_summary(cases):
	data = {}
	altered_cases = 0
	total_cases = 0
	total_valid_cases = 0
	data['altered_case_ids'] = []

	for case in cases:
		total_cases += 1
		# if data was not available, move to the next data point 
		#   without incrementing the count of valid data points
		if cases[case] == 'NA':
			continue
		else:
			val = int(cases[case])
			# increment the altered count by one for any valid non-zero result
			if int(val) not in [0, 1, -1]:
				altered_cases += 1
				data['altered_case_ids'].append(case)
			total_valid_cases += 1
	data['total_cases'] = total_cases
	data['total_valid_cases'] = total_valid_cases
	data['percent_decimal'] = float(altered_cases)/float(total_valid_cases)
	data['percent_display'] =  str( int( round(data['percent_decimal']*100, 0) ) ) + r'%'
	return data

def get_gene_mutation_data_summary(cases):
	'''Returns summary data gene mutation '''
	data = {}
	altered_cases = 0
	total_cases = 0
	total_valid_cases = 0
	data['altered_case_ids'] = []

	for case in cases:
		total_cases += 1
		total_valid_cases += 1
		if cases[case] == 'NaN' or str(cases[case]) == '0':
			continue
		else:
			data['altered_case_ids'].append(case)
			altered_cases += 1
	data['total_cases'] = total_cases
	data['total_valid_cases'] = total_valid_cases
	data['percent_decimal'] = float(altered_cases)/float(total_valid_cases)
	data['percent_display'] =  str( int( round(data['percent_decimal']*100, 0) ) ) + r'%'
	return data

def get_profile_data_API_response(url):
	'''This function returns an object that with same interface as an open file.'''
	response = requests.get(url)
	file_obj = StringIO.StringIO(response.text)
	return file_obj

def get_profile_data_file_name(genetic_profile_id, gene_list):
	'''This function returns the file_path to be used to locate a local version of the data returned by the API.
	It follows the file naming convention of simply using the URL-encoded querystring of the apply call (with .txt as an extension)'''
	file_dir = '../data_for_cbioportal_cmd_util/'
	payload = get_request_payload(genetic_profile_id, gene_list)
	query = 'cmd={cmd}&genetic_profile_id={genetic_profile_id}&id_type={id_type}&gene_list={gene_list}&case_set_id={case_set_id}'.format(**payload)
	file_name = file_dir + urllib.quote(query) + '.txt'
	return file_name

def get_profile_data_from_file(file_name):
	file_obj = open(file_name, 'r')
	return file_obj

def get_profile_data_URL(genetic_profile_id, gene_list):
	payload = get_request_payload(genetic_profile_id, gene_list)	
	return 'http://www.cbioportal.org/webservice.do?cmd={cmd}&genetic_profile_id={genetic_profile_id}&id_type={id_type}&gene_list={gene_list}&case_set_id={case_set_id}'.format(**payload)

def main():
	if ( len(sys.argv) == 1 ) or ( len(sys.argv) > 4 ):
		print('usage: gbm_summarize.py gene1 [gene2, gene3]')
	else:
		# get all arguments except the argv[0] (the script file), add them to a comma-delimited string
		gene_list = sys.argv[1:]
		display_gene_alterations_summary_descriptions(gene_list)

if __name__ == '__main__':
	main()