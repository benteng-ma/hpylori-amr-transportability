import csv
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
UPDATE=ROOT/'release_updates/v1.0.1'

def test_derived_audit_preserves_denominators_and_discrepancy_classes():
    with (UPDATE/'S12_provenance_audit.csv').open(encoding='utf-8',newline='') as f:
        rows=list(csv.reader(f))
    rows=rows[1:]
    assert len(rows)==57
    assert len({r[0] for r in rows})==57
    assert sum(r[21]=='yes' for r in rows)==57
    assert sum(r[22]=='yes' for r in rows)==31
    explained=[r for r in rows if r[23]=='Rule-definition explained: source rule includes N87Y; frozen panel excludes N87Y']
    assert len(explained)==1 and explained[0][0]=='SHZY43'
    assert sum(r[22]=='no' for r in rows)-len(explained)==25
    assert sum(r[24]=='FN' for r in rows)==16

def test_first_substantive_citations_are_sequential():
    text=(UPDATE/'manuscript_scientific_text.md').read_text(encoding='utf-8')
    for kind,count in [('Figure',12),('Table',20)]:
        spans=re.findall(rf'(?:Supplementary )?{kind}s?\s+(S\d+[a-d]?(?:\s*(?:[-–,]|and)\s*(?:S?\d+[a-d]?|[a-d]))*)',text)
        ordered=[]
        for span in spans:
            for n in re.findall(r'S(\d+)',span):
                n=int(n)
                if n not in ordered:ordered.append(n)
        assert ordered==list(range(1,count+1)),ordered

def test_number_map_is_complete_and_audit_is_not_new_validation():
    mapping=json.loads((UPDATE/'numbering_map.json').read_text(encoding='utf-8'))
    assert sorted(mapping['figures'].values())==list(range(1,13))
    assert sorted(mapping['tables'].values())==list(range(1,21))
    assert mapping['tables']['21']==12
    text=(UPDATE/'README.md').read_text(encoding='utf-8')
    assert 'does not add a cohort' in text
    assert 'not a new independent validation' in text
