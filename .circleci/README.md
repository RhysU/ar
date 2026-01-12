# CircleCI Configuration

This directory contains the CircleCI configuration for automated testing and
documentation deployment.

## Jobs

### `test`

Runs the test suite across multiple Python versions and NumPy versions.

### `deploy-docs`

Builds Doxygen documentation then deploys it to the `gh-pages` branch on GitHub
Pages.  This job only runs on the `master` branch after all tests pass.

## GitHub Pages Deployment Setup

To enable automatic documentation deployment, configure the following in your
CircleCI project settings:

### Required Environment Variables

Add these environment variables in CircleCI Project Settings → Environment
Variables:

1. **GH_NAME**: Your GitHub username or bot name for commits
   - Example: `CircleCI Documentation Bot`

2. **GH_EMAIL**: Email address for documentation commits
   - Example: `ci-bot@example.com`

3. **GH_TOKEN**: GitHub Personal Access Token with `repo` scope
   - Allows CircleCI to push to the `gh-pages` branch
   - Create at: https://github.com/settings/tokens
   - Required permissions: This repository only, Content, Read-Write

### GitHub Repository Setup

1. The `gh-pages` branch will be automatically created on first deployment

2. After the first deployment, enable GitHub Pages in repository settings:
   - Go to Settings → Pages
   - Source: Deploy from a branch
   - Branch: `gh-pages` / `root`

3. Your documentation will be available at:
   `https://<username>.github.io/<repository>/`

### Security Note

The `GH_TOKEN` should be kept secret. Never commit it to the repository.
CircleCI environment variables are encrypted and not exposed in build logs.

## References

- [CircleCI Environment Variables](https://circleci.com/docs/env-vars/)
- [GitHub Personal Access Tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
- [GitHub Pages Documentation](https://docs.github.com/en/pages)
