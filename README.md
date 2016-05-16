# cBioPortal Command Line Utilities

## Description
cBioPortal Command Line Utilities is project focused on providing simple command line tools for retrieving cancer genomic data using a web API and summarizing the results.

## Intstallation
### Python Package Dependencies
* [requests](https://pypi.python.org/pypi/requests/)

To install this package run the following command:
```bash
$ pip install requests
```
If that fails you may need to refer to the PIP [documentation](https://pip.pypa.io/en/stable/installing/).

## Background
This project utilize the cBioPortal web API. For more information on the project see the [documentation](http://www.cbioportal.org/web_api.jsp).

## Getting started

### Example: getting summary results for gene TP53

To run this analysis the following command can be used:
```bash
$ python gbm_summarize.py TP53
```

The command should return results in the following format:
```bash
TP53 is copy number altered in 23% of all cases.
TP53 is mutated in 29% of all cases.
```

### Example: getting summary results for multiple genes: TP53, MDM4, MDM2 

To run this analysis the following command can be used:
```bash
$ python gbm_summarize.py TP53 MDM2 MDM4
```

The command should return results in the following format:
```bash
Total % of cases where MDM4 is altered by either mutation or copy number alteration: 50% of all cases.
Total % of cases where MDM2 is altered by either mutation or copy number alteration: 40% of all cases.
Total % of cases where TP53 is altered by either mutation or copy number alteration: 30% of all cases.
```