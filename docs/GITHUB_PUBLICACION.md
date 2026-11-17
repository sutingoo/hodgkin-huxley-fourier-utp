# Publicación en GitHub (Ubuntu)

Este proyecto no se ha publicado automáticamente; el usuario deberá crear el repositorio, autenticar GitHub CLI y realizar `push`. Nombre sugerido: `sutingoo/hodgkin-huxley-fourier-utp`.

## Instalación de Git y GitHub CLI
```bash
sudo apt update
sudo apt install -y git gh
gh auth login
```

## Publicar desde la carpeta del proyecto
```bash
cd ~/Downloads/HHrepo
git init
git add .
git commit -m 'Reproducible HH Fourier study with Allen and OpenCOR evidence'
git branch -M main
gh repo create sutingoo/hodgkin-huxley-fourier-utp --public --source=. --remote=origin --push
```
Si el repositorio ya existe, usar `git remote add origin https://github.com/sutingoo/hodgkin-huxley-fourier-utp.git` y `git push -u origin main`, sin crear uno nuevo. **Comprueba antes los archivos versionados**: `git status`, `git ls-files`, `git ls-files | grep -E '\.(nwb|env)$'`.

## Después de publicar
1. Abrir el enlace real y comprobar que `README.md`, `src/`, `docs/`, `tests/` y `resultados/` sean visibles.
2. Reemplazar el texto `PENDIENTE DE PUBLICACIÓN` del Word por la URL real.
3. Guardar SHA del commit: `git rev-parse HEAD`; citarlo en la sección de reproducibilidad del Word.
4. No exponer tokens, contraseñas, archivos NWB enormes ni archivos de notas/ID estudiantiles innecesarios.
