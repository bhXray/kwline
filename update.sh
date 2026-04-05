#!/usr/bin/env bash
set -euo pipefail

REMOTE_NAME="origin"
REMOTE_URL="https://github.com/bhXray/kwline.git"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: this script must be run inside a Git repository."
  exit 1
fi

current_branch="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$current_branch" == "HEAD" ]]; then
  echo "Error: detached HEAD detected. Please checkout a branch first."
  exit 1
fi

if git remote get-url "$REMOTE_NAME" >/dev/null 2>&1; then
  existing_url="$(git remote get-url "$REMOTE_NAME")"
  if [[ "$existing_url" != "$REMOTE_URL" ]]; then
    echo "Updating $REMOTE_NAME URL to: $REMOTE_URL"
    git remote set-url "$REMOTE_NAME" "$REMOTE_URL"
  fi
else
  echo "Adding remote $REMOTE_NAME -> $REMOTE_URL"
  git remote add "$REMOTE_NAME" "$REMOTE_URL"
fi

commit_msg="${*:-chore: sync local changes $(date '+%Y-%m-%d %H:%M:%S')}"

git add -A

if git diff --cached --quiet; then
  echo "No staged changes to commit."
else
  git commit -m "$commit_msg"
fi

echo "Pushing branch '$current_branch' to $REMOTE_NAME ..."
git push -u "$REMOTE_NAME" "$current_branch"

echo "Done."