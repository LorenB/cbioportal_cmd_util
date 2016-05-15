import requests
import csv
import StringIO
import urllib
import os
import sys
import pprint



def get_case_set_data(gene_list):
	case_set_data = {'cases': {}, 'gene_meta_data':{}}
	genetic_profile_ids = ['gbm_tcga_gistic', 'gbm_tcga_mutations']
	for genetic_profile_id in genetic_profile_ids:
		add_genetic_profile(case_set_data, genetic_profile_id, gene_list)
		add_genetic_profile(case_set_data, genetic_profile_id, gene_list)	
	return case_set_data

def display_case_set_summary(case_set):
	display = [{'analysis': 'profile_gbm_tcga_mutations_summary', 'text': '{gene} is mutated in {result_text}% of all cases.', 'summary_function': get_mutation_percent},
				{'analysis': 'profile_gbm_tcga_gistic_summary',  'text': '{gene} is copy number altered in {result_text}% of all cases.', 'summary_function': get_copy_number_altered_percent},
				{'analysis': 'aggregate', 'text': 'Total % of cases where {gene} is altered by either mutation or copy number alteration: {result_text}% of all cases.', 'summary_function': get_all_alterated_percent}]
	
	case_set_summary = get_case_set_summary(case_set)
	tokens = {}
	for gene in case_set_summary:
		tokens = {'gene': gene}
		for metric in display:
			result_percent= metric['summary_function'](case_set_summary, gene)
			tokens['result_text'] = '{0:.{1}f}'.format(result_percent*100, 0)
			print metric['text'].format(**tokens)

def get_mutation_percent(case_set_summary, gene):
	return float(case_set_summary[gene]['mutated_case_count']) / float(case_set_summary[gene]['total_case_count'])

def get_copy_number_altered_percent(case_set_summary, gene):
	return float(case_set_summary[gene]['copy_number_alterated_case_count']) / float(case_set_summary[gene]['total_case_count'])

def get_all_alterated_percent(case_set_summary, gene):
	return float(case_set_summary[gene]['copy_number_alterated_case_count'] + case_set_summary[gene]['mutated_case_count'] - case_set_summary[gene]['multiple_alterations_case_count']) / float(case_set_summary[gene]['total_case_count'])

def is_mutated(case, gene):
	genetic_profile_id = 'gbm_tcga_mutations'
	profile_result = case[gene][genetic_profile_id] 
	if profile_result == 'NaN' or str(profile_result) == '0':
		return False
	else:
		return True

def is_copy_number_alterated(case, gene):
	genetic_profile_id = 'gbm_tcga_gistic'
	profile_result = case[gene][genetic_profile_id]
	# @TODO: determine best way to handle 'NA' values
	if int(profile_result) not in [0, 1, -1]:
		return True
	else:
		return False

def get_case_set_summary(case_set):
	case_set_summary = { }
	total_case_count = len(case_set['cases'])
	for gene in case_set['gene_meta_data']:
		mutated_case_count = 0
		copy_number_alterated_case_count = 0
		multiple_alterations_case_count = 0

		for case in case_set['cases']:
			if is_copy_number_alterated(case_set['cases'][case], gene):
				copy_number_alterated_case_count += 1
			if is_mutated(case_set['cases'][case], gene):
				mutated_case_count += 1
			if is_copy_number_alterated(case_set['cases'][case], gene) and is_mutated(case_set['cases'][case], gene):
				multiple_alterations_case_count += 1
		
		mutated_case_percent = float(mutated_case_count) / float(total_case_count)
		copy_number_alterated_percent = float(copy_number_alterated_case_count) / float(total_case_count)
		all_alterated_percent = float(copy_number_alterated_case_count + mutated_case_count - multiple_alterations_case_count) / float(total_case_count)
		case_set_summary[gene] = {'total_case_count': total_case_count,
									'mutated_case_count': mutated_case_count,
									'copy_number_alterated_case_count': copy_number_alterated_case_count,
									'multiple_alterations_case_count': multiple_alterations_case_count,
									'all_alterated_percent': all_alterated_percent}

	return case_set_summary

def add_genetic_profile(case_set_data, genetic_profile_id, gene_list, useLocalFile=False):
	'''Returns one genetic profile for ever gene provided.'''
	gene_alterations = {'meta_data': [], 'gene_meta_data': {}, 'cases': {} }
	gene_meta_data_fields = ['COMMON', 'GENE_ID']
	
	if type(gene_list) == list():
		gene_list = ','.join(gene_list)
	
	if useLocalFile:
		file_obj = get_profile_data_from_file(get_profile_data_file_name(genetic_profile_id, gene_list))
	else:
		url = get_profile_data_URL(genetic_profile_id, gene_list)
		# create a object with same interface as open file to allow for the use of file-based parsing methods
		file_obj = get_profile_data_API_response(url)

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
				if not (gene in case_set_data['gene_meta_data']):
					case_set_data['gene_meta_data'][gene] = {}
				case_set_data['gene_meta_data'][gene][field] = row[field]
			else:
				if not (field in case_set_data['cases']):
					case_set_data['cases'][field] = {}
				if not (gene in case_set_data['cases'][field]):
					case_set_data['cases'][field][gene] = {}
				case_set_data['cases'][field][gene][genetic_profile_id] = row[field]
	file_obj.close()


def get_request_payload(genetic_profile_id, gene_list):
	payload = {}
	payload['cmd'] = 'getProfileData'
	payload['id_type'] = 'gene_symbol'
	payload['case_set_id'] = 'gbm_tcga_cnaseq'
	payload['genetic_profile_id'] = genetic_profile_id
	payload['gene_list'] = gene_list
	return payload

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
		gene_list = ','.join(sys.argv[1:])
		case_set_data = get_case_set_data(gene_list)
		display_case_set_summary(case_set_data)

if __name__ == '__main__':
	main()