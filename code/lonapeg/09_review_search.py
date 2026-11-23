"""LONAPEG-FAERS Part B: literature searches (PubMed E-utilities, Europe PMC incl. preprints/abstracts). Search date recorded."""
import json, urllib.request, urllib.parse, time, datetime, os
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','review','search')
PUBMED='(lonapegsomatropin[tiab] OR Skytrofa[tiab] OR "TransCon GH"[tiab] OR "TransCon hGH"[tiab] OR "TransCon growth hormone"[tiab] OR "ACP-011"[tiab] OR "ACP-001"[tiab] OR lonapegsomatropin[nm])'
EPMC='(lonapegsomatropin OR skytrofa OR "TransCon GH" OR "TransCon hGH" OR "TransCon growth hormone" OR "ACP-011")'
def get(u): return json.load(urllib.request.urlopen(u, timeout=60))
r=get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=1000&retmode=json&term="+urllib.parse.quote(PUBMED))['esearchresult']
ids=r['idlist']; recs=[]
for i in range(0,len(ids),150):
    s=get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="+",".join(ids[i:i+150]))['result']
    for k in ids[i:i+150]:
        x=s[k]; recs.append(dict(source='PubMed',pmid=k,title=x['title'],journal=x['source'],year=x['pubdate'][:4],pubtype='; '.join(x.get('pubtype',[])),
            doi=next((a['value'] for a in x.get('articleids',[]) if a['idtype']=='doi'),''),pmcid=next((a['value'] for a in x.get('articleids',[]) if a['idtype']=='pmc'),''),
            first_author=(x['authors'][0]['name'] if x.get('authors') else '')))
    time.sleep(0.4)
ep=[]; cursor='*'
while True:
    d=get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=1000&resultType=lite&cursorMark="+urllib.parse.quote(cursor)+"&query="+urllib.parse.quote(EPMC))
    for x in d['resultList']['result']:
        ep.append(dict(source='EuropePMC:'+x.get('source',''),pmid=x.get('pmid',''),title=x.get('title',''),journal=x.get('journalTitle',''),year=x.get('pubYear',''),
                       pubtype=x.get('pubType',''),doi=x.get('doi',''),pmcid=x.get('pmcid',''),first_author=(x.get('authorString','').split(',')[0])))
    nc=d.get('nextCursorMark')
    if not nc or nc==cursor or not d['resultList']['result']: break
    cursor=nc
meta=dict(date=datetime.datetime.now(datetime.timezone.utc).isoformat(), pubmed_query=PUBMED, pubmed_count=int(r['count']), epmc_query=EPMC, epmc_count=len(ep))
json.dump(dict(meta=meta,pubmed=recs,epmc=ep),open(f'{OUT}/search_results.json','w'),indent=1); print(meta)
