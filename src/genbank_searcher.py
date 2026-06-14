import os

from Bio import Entrez
import os

Entrez.email = "bcominscheffel@gmail.com"

# Query para buscar genes TEF1-alpha em fungos com cds completa ou parcial diretamente do banco de dados GenBank, utilizando a API do NCBI.

# CDS (complete coding sequence) é a sequência de DNA que codifica uma proteína completa, incluindo os códons de início e término.

QUERY = (
    'HIV-1[Organism] AND (pol[Gene] OR RT[Gene] OR protease[Gene]) '
    'AND "drug resistance"[Title/Abstract] '
    'AND "population sequencing"[All Fields] OR "Sanger"[All Fields] '
    'AND 500:3500[Sequence Length] '
    'NOT RefSeq[Keyword] NOT clone[Title]'
)

MAX_RECORDS = 1000
BATCH_SIZE = 50

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

OUTPUT_FILE = os.path.join(BASE_DIR, "../assets/genbank_data/hiv_env_teste.gb")

def main():
    print("Buscando genes HIV env...")

    handle = Entrez.esearch(
        db="nucleotide",
        term=QUERY,
        retmax=MAX_RECORDS
    )
    result = Entrez.read(handle)
    ids = result["IdList"] #type: ignore
    handle.close()

    print(f"{len(ids)} IDs selecionados")

    with open(OUTPUT_FILE, "w", encoding="utf-8") as out:
        for i in range(0, len(ids), BATCH_SIZE):
            batch = ids[i:i + BATCH_SIZE]

            fetch_handle = Entrez.efetch(
                db="nucleotide",
                id=batch,
                rettype="gb",
                retmode="text"
            )

            records = fetch_handle.read()  # lê o texto do GenBank
            fetch_handle.close()
            out.write(records)

            print(f"{min(i + BATCH_SIZE, len(ids))}/{len(ids)} registros baixados")

    print("Arquivo final pronto")

if __name__ == "__main__":
    main()
