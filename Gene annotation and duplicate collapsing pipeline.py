#!/usr/bin/env python3
"""Gene annotation and duplicate collapsing pipeline."""

import argparse
import gzip
import json
import logging
import re
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import mygene
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

ENSEMBL_REST = "https://rest.ensembl.org"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INPUT_FILENAME = "top_genes_mapped.tsv"
ANNOTATED_SUFFIX = "_annotated.tsv"
COLLAPSED_SUFFIX = "_annotated_collapsed.tsv"
FAILED_MAPPINGS_FILE = "failed_mappings.tsv"
DEFAULT_CACHE_NAME = "gene_annotation_cache.sqlite3"
NCBI_TOOL = "gene-annotation-pipeline"
DEFAULT_NCBI_EMAIL = ""

def utc_now():
    return datetime.now(timezone.utc).isoformat()

def make_session():
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.headers.update({
        "Content-Type": "application/json",
        "User-Agent": "gene-annotation-pipeline/1.0",
    })
    return session

SESSION = make_session()

    }
    )