import requests
import time
import sys # Pour afficher les erreurs sur stderr

# --- Configuration ---
input_filename = "data/GWAS_sign_SNP.txt"
output_filename = "data/snp_coordinates_grch37.tsv"
api_url_base = "http://grch37.rest.ensembl.org/variation/homo_sapiens"
# Temps de pause entre les requêtes (en secondes) pour respecter le serveur
sleep_time = 0.5
# --- Fin Configuration ---

# Compteurs pour le résumé
count_success = 0
count_fail = 0

# Essayer d'ouvrir les fichiers
try:
    infile = open(input_filename, 'r')
except FileNotFoundError:
    print(f"ERREUR: Le fichier d'input '{input_filename}' n'a pas été trouvé!", file=sys.stderr)
    sys.exit(1) # Arrêter le script

try:
    outfile = open(output_filename, 'w')
except IOError:
    print(f"ERREUR: Impossible d'écrire dans le fichier de sortie '{output_filename}'!", file=sys.stderr)
    infile.close()
    sys.exit(1)


# Écrire l'en-tête dans le fichier de sortie
outfile.write("rsID\tChromosome\tPosition_GRCh37\n")

print(f"Lecture des rsIDs depuis '{input_filename}'...")
print(f"Écriture des coordonnées GRCh37 dans '{output_filename}'...")

# Lire chaque ligne (rsID) du fichier d'input
for line in infile:
    rsid = line.strip() # Enlever les espaces/retours à la ligne

    if not rsid: # Ignorer les lignes vides
        continue

    print(f"Traitement de {rsid} ... ", end="")

    request_url = f"{api_url_base}/{rsid}?content-type=application/json"

    chromosome = "NA"
    position = "NA"

    try:
        # Faire la requête à l'API Ensembl
        response = requests.get(request_url, headers={"Content-Type": "application/json"}, timeout=10) # Timeout de 10s

        # Vérifier si la requête a échoué (ex: 404 Not Found si rsid inconnu)
        response.raise_for_status() # Lève une exception pour les erreurs HTTP

        # Si la requête réussit (status code 200)
        data = response.json()

        # Essayer d'extraire les coordonnées GRCh37
        # On essaie d'abord au niveau racine, puis dans 'mappings' si nécessaire
        try:
            chromosome = data['seq_region_name']
            position = data['start']
        except KeyError:
             # Si pas au niveau racine, essayer dans le premier mapping
            try:
                mapping = data.get('mappings') # Utilise .get pour éviter KeyError si 'mappings' n'existe pas
                if mapping and isinstance(mapping, list) and len(mapping) > 0:
                     chromosome = mapping[0]['seq_region_name']
                     position = mapping[0]['start']
                else:
                     # Si 'mappings' est absent ou vide, on ne trouve pas les coordonnées
                      raise KeyError # Provoque le passage au bloc except extérieur
            except (KeyError, IndexError, TypeError):
                print("Échec (Structure JSON inattendue ou coordonnées manquantes).")
                count_fail += 1
                outfile.write(f"{rsid}\tPARSE_ERROR\tPARSE_ERROR\n")
                time.sleep(sleep_time) # Pause même en cas d'erreur de parsing
                continue # Passer au rsID suivant

        # Si on a réussi à extraire chromosome et position
        print(f"OK ({chromosome}:{position})")
        outfile.write(f"{rsid}\t{chromosome}\t{position}\n")
        count_success += 1

    except requests.exceptions.HTTPError as http_err:
        # Gérer les erreurs HTTP (ex: 404 Not Found)
        print(f"Échec (Erreur HTTP: {http_err})")
        count_fail += 1
        # Pas besoin d'écrire ici, car 'NA' sont les valeurs par défaut

    except requests.exceptions.RequestException as req_err:
        # Gérer les autres erreurs de requête (connexion, timeout...)
        print(f"Échec (Erreur Requête: {req_err})")
        count_fail += 1
        # Pas besoin d'écrire ici

    # Écrire les NA si une erreur HTTP ou RequestException s'est produite avant l'extraction
    if chromosome == "NA":
         outfile.write(f"{rsid}\t{chromosome}\t{position}\n")

    # Pause pour être respectueux envers le serveur Ensembl
    time.sleep(sleep_time)

# Fermer les fichiers
infile.close()
outfile.close()

print("\n-----------------------------------------------------")
print("Terminé.")
print(f"Succès: {count_success}")
print(f"Échecs/Non trouvés/Erreurs: {count_fail}")
print(f"Résultats sauvegardés dans : '{output_filename}'")
print("Veuillez vérifier le fichier de sortie.")
print("-----------------------------------------------------")