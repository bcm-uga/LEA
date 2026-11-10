# Revue du correctif de sécurité ANCESTRYMAP

Date : 10 novembre 2026.
Version : LEA 3.25.0.
Branche de travail : `fix/ancestrymap-buffer-overflow`.
Référence de comparaison : `devel`.

## Périmètre

Revue des modifications non commitées de `src/convert/ancestrymap.c` et `src/convert/ancestrymap.h` par rapport à `devel` : correction des débordements, bornes des buffers, traitement des lectures et de la fin de fichier, comptage des individus, indexation de la matrice et conservation du comportement des conversions valides.

La revue repose sur l'inspection du diff et du code, y compris les fonctions d'entrée/sortie appelées. Aucun nouveau test n'a été exécuté pendant cette revue. Les résultats des tests précédemment exécutés sont documentés dans `ANCESTRYMAP_SECURITY_FIX_REPORT.md`.

## Constat actionnable

### Sévérité faible — Fuite de descripteur lors du rejet d'un nom de SNP

**Caractère : non bloquant.**

**Références :** `src/convert/ancestrymap.c:227` ; erreurs levées par le validateur aux lignes 37 et 40.

Dans `nb_ind_ancestrymap`, le fichier est ouvert à la ligne 222, puis `check_snp_name` est appelé à la ligne 227. Si le premier champ est absent ou contient au moins 512 octets, le validateur déclenche une erreur R. Ce transfert de contrôle évite le `fclose(File)` situé à la ligne 238.

Des appels répétés avec des fichiers invalides peuvent donc accumuler des descripteurs ouverts et, à terme, empêcher l'ouverture de fichiers pour des conversions valides. Ce constat découle du chemin de contrôle ; l'épuisement des descripteurs n'a pas été reproduit pendant la revue.

**Correction recommandée :** fermer le flux de la passe de comptage avant de lever une erreur de validation, ou prévoir un mécanisme de nettoyage exécuté lors d'une erreur R. Cette amélioration peut faire l'objet d'un suivi distinct.

## Décision

Aucun problème bloquant de correction, de bornes des buffers, de gestion de la fin de fichier, de comptage des individus ou d'indexation de la matrice n'a été identifié dans les modifications examinées. Aucune régression du comportement des conversions ANCESTRYMAP valides n'a été identifiée.

**Le correctif est prêt à être commité**, avec le suivi non bloquant décrit ci-dessus.

Aucun fichier source n'a été modifié pour cette revue ou pour la rédaction de ce rapport. Aucun commit ni push n'a été effectué.
