#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Edmund Su
    Linda Xiang
"""

import os
import sys
import argparse
import subprocess
import json
import re
import hashlib
import uuid
import tarfile
from datetime import date
import copy
from glob import glob
import yaml
import csv
import io
import shutil
#LUCA-KR.DO231106.SA602282.wxs.20210112.gatk-mutect2.somatic.snv.open-filter.vcf.gz
workflow_process_map = {
    'DNA Seq Germline Workflow': 'snv'
}

tool_list = ['strelka']

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()

def get_files_info(file_to_upload, date_str, analysis_dict, process_indicator,tool,new_dir):
    file_info = {
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'info': {
            'data_category': "Simple Nucleotide Variation",
        }
    }
    ### deepvariant
    if tool=="deepvariant":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['DeepVariant']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['DeepVariant']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    elif tool=="strelka":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['Strelka']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['Strelka']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    elif tool=="tiddit":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['Tiddit']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['Tiddit']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    elif tool=="haplotypecaller":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['haplotypecaller']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['haplotypecaller']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    elif tool=="manta":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['Manta']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['Manta']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    elif tool=="freebayes":
        if re.match(r'.*.vcf.gz$', file_to_upload):
            file_type = 'VCF'
            file_info.update({'dataType': 'Raw SNV Calls'})
            file_info['info'].update({'analysis_tools': ['Freebayes']})
        elif re.match(r'.*.vcf.gz.tbi$', file_to_upload):
            file_type = 'TBI'
            file_info.update({'dataType': 'VCF Index'})
            file_info['info'].update({'analysis_tools': ['Freebayes']})
        else:
            sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)
    #elif tool=="cnvkit":
    else:
        sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)

    #LUCA-KR.DO231106.SA602282.wxs.20210112.gatk-mutect2.somatic.snv.open-filter.vcf.gz.tbi"
    #"TEST-PR.DO250183.SA610228.wxs.20230501.snv-strelka.gvcf.gz",
    suffix={
        "VCF":"vcf.gz",
        "TBI": "vcf.gz.tbi",
    }
    # file naming patterns:
    #   pattern:  <argo_study_id>.<argo_donor_id>.<argo_sample_id>.<experiment_strategy>.<date>.<process_indicator>.<file_type>.<file_ext>
    #   process_indicator: pre-alignment, alignment(aligner), post-alignment(caller)
    #   example: TEST-PR.DO250183.SA610229.rna-seq.20200319.star.genome_aln.cram
    new_fname = '.'.join([
        analysis_dict['studyId'],
        analysis_dict['samples'][0]['donor']['donorId'],
        analysis_dict['samples'][0]['sampleId'],
        analysis_dict['experiment']['experimental_strategy'].lower() if analysis_dict['experiment'].get('experimental_strategy') else analysis_dict['experiment']['library_strategy'],
        date_str,
        process_indicator,
        suffix[file_type]
      ])    
    
    file_info['fileName'] = new_fname
    file_info['fileType'] = file_type

    if re.match(r'cnvkit', file_to_upload):
        with tarfile.open(file_to_upload, 'r') as tar:
            for member in tar.getmembers():
                file_info['info']['files_in_tgz'].append(member.name)

            new_dir = 'out'
            try:
                os.mkdir(new_dir)
            except FileExistsError:
                pass

            dst = os.path.join(os.getcwd(), new_dir, new_fname)
            os.symlink(os.path.abspath(file_to_upload), dst)
    else:
        shutil.copyfile(os.path.realpath(file_to_upload),"/".join([new_dir,new_fname]))
        ##os.symlink(os.path.realpath(file_to_upload),"/".join([new_dir,new_fname]))

    return file_info

def get_basename(metadata):
    study_id = metadata['studyId']
    donor_id = metadata['samples'][0]['donor']['donorId']
    sample_id = metadata['samples'][0]['sampleId']

    if not sample_id or not donor_id or not study_id:
      sys.exit('Error: missing study/donor/sample ID in the provided metadata')

    return ".".join([study_id, donor_id, sample_id])

def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
      for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
        sample.pop(item, None)
        sample['specimen'].pop(item, None)
        sample['donor'].pop(item, None)

    return samples

def prepare_tarball(sampleId, qc_files, tool_list):

    tgz_dir = 'tarball'
    try:
      os.mkdir(tgz_dir)
    except FileExistsError:
      pass

    files_to_tar = {}
    for tool in tool_list:
      if not tool in files_to_tar: files_to_tar[tool] = []
      for f in sorted(qc_files):
        if tool in f:
          files_to_tar[tool].append(f)   
    
    for tool in tool_list:
      if not files_to_tar[tool]: continue
      tarfile_name = f"{tgz_dir}/{sampleId}.{tool}.tgz"
      with tarfile.open(tarfile_name, "w:gz", dereference=True) as tar:
        for f in files_to_tar[tool]:
          tar.add(f, arcname=os.path.basename(f))

def main():
    """
    Python implementation of tool: payload-gen-qc
    """

    parser = argparse.ArgumentParser(description='Tool: payload-gen-qc')
    parser.add_argument("-a", "--metatada-analysis", dest="metadata_analysis", required=True,
                        help="Input metadata analysis", type=str)
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="All files to upload")
    parser.add_argument("-g", "--genome_annotation", dest="genome_annotation", default="", help="Genome annotation")
    parser.add_argument("-b", "--genome_build", dest="genome_build", default="", help="Genome build")
    parser.add_argument("-w", "--wf-name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-s", "--wf-session", dest="wf_session", required=True, help="workflow session ID")
    parser.add_argument("-v", "--wf-version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-p", "--pipeline_yml", dest="pipeline_yml", required=False, help="Pipeline info in yaml")
    parser.add_argument("-t", "--tool", dest="tool", required=True,type=str, help="Tool used for variant calling",
    choices=['strelka','cnvkit','deepvariant','tiddit','manta','haplotypecaller','freebayes']
    )

    args = parser.parse_args()
    
    with open(args.metadata_analysis, 'r') as f:
      analysis_dict = json.load(f)

    pipeline_info = {}
    if args.pipeline_yml:
      with open(args.pipeline_yml, 'r') as f:
        pipeline_info = yaml.safe_load(f)

    payload = {
        'analysisType': {
            'name': 'variant_processing'
        },
        'studyId': analysis_dict.get('studyId'),
        'info': {},
        'workflow': {
            'workflow_name': args.wf_name,
            'workflow_version': args.wf_version,
            'session_id': args.wf_session,
            'inputs': [
                {
                    'analysis_type': analysis_dict['analysisType']['name'],
                    'input_analysis_id': analysis_dict.get('analysisId')
                }
            ],
            'info': pipeline_info
        },
        'files': [],
        'experiment': analysis_dict.get('experiment'),
        'samples': get_sample_info(analysis_dict.get('samples'))
    }
    if args.genome_build:
      payload['workflow']['genome_build'] = args.genome_build
    if args.genome_annotation:
      payload['workflow']['genome_annotation'] = args.genome_annotation

    # pass `info` dict from seq_experiment payload to new payload
    if 'info' in analysis_dict and isinstance(analysis_dict['info'], dict):
      payload['info'] = analysis_dict['info']
    else:
      payload.pop('info')

    if 'library_strategy' in payload['experiment']:
      experimental_strategy = payload['experiment'].pop('library_strategy')
      payload['experiment']['experimental_strategy'] = experimental_strategy

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    # generate date string
    date_str = date.today().strftime("%Y%m%d")

    # prepare tarball to include all QC files generated by one tool
    ##prepare_tarball(analysis_dict['samples'][0]['sampleId'], args.files_to_upload, tool_list)

    process_indicator = ".".join([args.tool,"germline",workflow_process_map.get(args.wf_name)])
    for f in sorted(args.files_to_upload):
      file_info = get_files_info(f, date_str, analysis_dict, process_indicator,args.tool,new_dir)
      payload['files'].append(file_info)

    with open("%s.%s.payload.json" % (str(uuid.uuid4()), args.wf_name.replace(" ","_")), 'w') as f:
        f.write(json.dumps(payload, indent=2))



if __name__ == "__main__":
    main()
