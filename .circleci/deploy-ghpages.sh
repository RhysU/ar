#!/bin/bash
# Deploy documentation to gh-pages branch
# Based on: https://github.com/DevProgress/onboarding/wiki/Using-Circle-CI-with-Github-Pages-for-Continuous-Delivery

# Exit with nonzero exit code if anything fails or is undefined
set -eu

# Get the directory containing the built documentation
SOURCE_DIR="${1:-docs/html}"
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

# Configure git
git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"

# Create a temporary directory and initialize a fresh git repository
TEMP_DIR=$(mktemp -d)
trap "rm -rf '$TEMP_DIR'" EXIT
cd "$TEMP_DIR"

# Initialize a new repository with gh-pages as the default branch
git init -b gh-pages

# Copy documentation into the repository
cp -R "${CIRCLE_WORKING_DIRECTORY}/${SOURCE_DIR}"/* .

# Create .nojekyll file to allow files/dirs with underscores
touch .nojekyll

# Commit all content
git add -A
git commit -m "Deploy documentation to GitHub Pages [ci skip]

Built from commit ${CIRCLE_SHA1} on branch ${CIRCLE_BRANCH}"

# Set remote and force push to replace gh-pages branch
git remote add origin "https://${GH_TOKEN}@github.com/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}.git"

echo "Pushing to gh-pages branch..."
git push --force --quiet origin gh-pages > /dev/null 2>&1

echo "Documentation deployed successfully to gh-pages"
