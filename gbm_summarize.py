import requests
import csv
import StringIO
import urllib
import os
import sys

def display_all_profile_summary_data(gene_list):
	genetic_profile_id_config = {
		'gbm_tcga_mutations': {'summary_function': get_mutated_percent_data, 'display_text': '{COMMON} is mutated in {value_display} of all cases.'},
		'gbm_tcga_gistic': {'summary_function': get_percent_of_cases_copy_number_altered_data, 'display_text': '{COMMON} is copy number altered in {value_display} of all cases.'}
	}
	summary_function = ''
	results = ''
	results_text = {}

	for genetic_profile_id in genetic_profile_id_config:
		summary_function = genetic_profile_id_config[genetic_profile_id]['summary_function']
		results =  get_profile_data(genetic_profile_id, gene_list, summary_function)
		for result in results:
			if result['COMMON'] not in results_text.keys():
				results_text[result['COMMON']] = []
			results_text[result['COMMON']].append( genetic_profile_id_config[genetic_profile_id]['display_text'].format(**result) )
	for gene_result_text in results_text:
		for result_line in results_text[gene_result_text]:
			print result_line

def get_profile_data(genetic_profile_id, gene_list, summary_function, useLocalFile=False):
	if isinstance(gene_list, list):
		gene_list_str = ','.join(gene_list_str)
	else:
		gene_list_str = gene_list
	url = get_profile_data_URL(genetic_profile_id, gene_list_str)
	
	if useLocalFile:
		file_obj = get_profile_data_from_file(get_profile_data_file_name(genetic_profile_id, gene_list))
	else:
		# create a object with same interface as open file to allow for the use of file-based parsing methods
		file_obj = get_profile_data_API_response(url)

	# create an array to store the position of each line the file, to allow for multiple rows of meta-data to proceed the column headers
	meta_data = []
	header_row_begin_position = 0
	for line in iter(file_obj.readline, ''):
		# check if the line begins with "#" as this is the convention for providing meta-data in the text returned by the API
		if line[0] == '#':
			meta_data.append(line)
			header_row_begin_position = file_obj.tell()
		else:
			break

	# get the position before the header row (the file curesor will be at the end of the end of the header row when the loop exits)
	file_obj.seek(header_row_begin_position)

	reader = csv.DictReader(file_obj, delimiter='\t')
	results = []
	for row in reader:
		results.append( summary_function(row) )
	file_obj.close()
	return results

def get_request_payload(genetic_profile_id, gene_list):
	payload = {}
	payload['cmd'] = 'getProfileData'
	payload['id_type'] = 'gene_symbol'
	payload['case_set_id'] = 'gbm_tcga_cnaseq'
	payload['genetic_profile_id'] = genetic_profile_id
	payload['gene_list'] = gene_list
	return payload

def get_percent_of_cases_copy_number_altered_data(profile):
	data = {}	
	altered_count = 0
	total_valid = 0

	meta_data_headers = ['COMMON', 'GENE_ID']

	for key in profile:
		# print(key, profile[key])
		if key in meta_data_headers:
			data[key] = profile[key]
			continue
		
		# if data was not available, move to the next data point 
		#   without incrementing the count of valid data points
		if profile[key] == 'NA':
			continue
		else:
			val = int(profile[key])
			# increment the altered count by one for any valid non-zero result
			if int(val) not in [0, 1, -1]:
				altered_count += 1
			total_valid += 1
	data['value'] = float(altered_count)/float(total_valid)
	data['value_display'] =  str( int( round(data['value']*100, 0) ) ) + r'%'	
	return data

def get_mutated_percent_data(profile):
	data = {}
	mutation_count = 0
	total_valid = 0
	meta_data_headers = ['COMMON', 'GENE_ID']
	for key in profile:
		if key in meta_data_headers:
			data[key] = profile[key]
			continue
		total_valid += 1
		if profile[key] == 'NaN' or str(profile[key]) == '0':
			continue
		else:
			mutation_count += 1
	data['value'] = float(mutation_count)/float(total_valid)
	data['value_display'] =  str( int( round(data['value']*100, 0) ) ) + r'%'
	return data	


def get_profile_data_API_response(url):
	# method returns an object that with same interface as an open file
	response = requests.get(url)
	file_obj = StringIO.StringIO(response.text)
	return file_obj



def get_profile_data_file_name(genetic_profile_id, gene_list):
	'''This function returns the file_path to be used to locate a local version of the data returned by the API.
	It follows the file naming convention of simply using the URL-encoded querystring of the apply call (with .txt as an extension)'''
	file_dir = './'
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
		display_all_profile_summary_data( gene_list )

if __name__ == '__main__':
	main()