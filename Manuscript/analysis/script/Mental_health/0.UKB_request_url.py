import pandas as pd
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import re


def parse_table_data(url):
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to get content from url")
        return None
    soup = BeautifulSoup(response.content, 'html.parser')
    tables = soup.find_all('table')
    rows = tables[1].find_all('tr')
    row0 = rows[0].find_all('td')
    row1 = rows[1].find_all('td')
    row2 = rows[2].find_all('td')
    res_lst = []
    res_lst.extend([row0[0].text.strip(), row0[1].text.strip()])
    res_lst.extend([row1[0].text.strip(), row1[1].text.strip()])
    res_lst.extend([row0[3].text.strip(), row0[4].text.strip()])
    res_lst.extend([row0[6].text.strip(), row0[7].text.strip()])
    res_lst.extend([row2[2].text.strip(), row2[3].text.strip()])
    divs = soup.find_all('div', class_='tabbertab')
    pattern = r"Data-Coding\s+[0-9]+"
    for div in divs:
        if div and div.find('h2', string='Data'):
            matches = re.findall(pattern, div.text)
            res_lst.extend(matches)
    return res_lst


df = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/trait.txt", sep="\t")
Field_lst = list(df.iloc[:, 0])
name_lst = list(df.iloc[:, 1])

with open("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/trait.url.txt", "w") as fout:
    for (field, name) in tqdm(zip(Field_lst, name_lst), total=len(Field_lst)):
        url = "https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={}".format(field)
        try:
            res_lst = parse_table_data(url)
            fout.write(f"{field}\t{name}\t" + "\t".join(res_lst) + "\n")
        except Exception as e:
            fout.write(f"{field}\t{name}\n")
