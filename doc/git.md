##### Publication après des modifications locales
1. Ajout des fichiers modifiés pour publication
```
git add file_name1 file_name2 ...
```
2. Commit
```
git commit -m "message explicatif"
```
3. Récupération d'éventuelles modifications faites par d'autres (peut nécessiter la résolution de conflits)
```
git pull
```
4. Publication de ses changements
```
git push
```

##### Forcer l'ajout de fichiers exclus par le .gitignore (pour les .rds par exemple)
```
git add -f file_name
```

##### Créer une nouvelle branche
```
git checkout -b branch_name
```

##### Changer de branche (existante)
```
git checkout branch_name
```

##### Oublier les modifications locales faites sur un fichier et revenir à la version distante
```
git checkout -- file_name
```

##### Importer le travail fait sur une branche dans la branche main
1. Se placer dans la branche main
```
git checkout main
```
2. Effectuer le merge (peut nécessiter la résolution de conflits)
```
git merge branch_name
```

##### Résoudre des conflits
1. Modifier les fichiers en conflits
2. S'assurer que la résolution des conflits fonctionne
3. Publier la résolution des conflits (voir ci-dessus)