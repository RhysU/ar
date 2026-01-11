#!/bin/bash
# Deploy documentation to gh-pages branch
# Based on: https://github.com/DevProgress/onboarding/wiki/Using-Circle-CI-with-Github-Pages-for-Continuous-Delivery

set -e # Exit with nonzero exit code if anything fails

# Get the directory containing the built documentation
SOURCE_DIR="${1:-docs/html}"

if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

# Configure git
git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"

# Clone the existing gh-pages branch into a temporary directory
TEMP_DIR=$(mktemp -d)
cd "$TEMP_DIR"

git clone --quiet --branch=gh-pages "https://${GH_TOKEN}@github.com/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}.git" gh-pages > /dev/null 2>&1 || {
    # If gh-pages doesn't exist, create it
    git clone --quiet "https://${GH_TOKEN}@github.com/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}.git" gh-pages > /dev/null 2>&1
    cd gh-pages
    git checkout --orphan gh-pages
    git rm -rf . > /dev/null 2>&1 || true
}

cd gh-pages

# Remove all existing content (except .git)
find . -maxdepth 1 ! -name '.git' ! -name '.' ! -name '..' -exec rm -rf {} \; 2>/dev/null || true

# Copy new documentation
cp -R "${CIRCLE_WORKING_DIRECTORY}/${SOURCE_DIR}"/* .

# Create .nojekyll file to allow files/dirs with underscores
touch .nojekyll

# Check if there are changes to commit
if git diff --quiet && git diff --cached --quiet; then
    echo "No changes to documentation"
    cd "$CIRCLE_WORKING_DIRECTORY"
    rm -rf "$TEMP_DIR"
    exit 0
fi

# Commit and push
git add -A
git commit -m "Deploy documentation to GitHub Pages [ci skip]

Built from commit ${CIRCLE_SHA1} on branch ${CIRCLE_BRANCH}"

echo "Pushing to gh-pages branch..."
git push --force --quiet origin gh-pages > /dev/null 2>&1

echo "Documentation deployed successfully to gh-pages"

# Cleanup
cd "$CIRCLE_WORKING_DIRECTORY"
rm -rf "$TEMP_DIR"
