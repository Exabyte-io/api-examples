#!/bin/bash

# Change branch name in the Colab URLs in .ipynb files

if [ -z "$1" -o -z "$2" ]; then
    echo "Usage: change-branch-in-urls.sh <old-branch-name> <new-branch-name>"
    exit 1
fi

old_branch_name="$1"
new_branch_name="$2"

_RED=$(tput setaf 1; tput setab 7)
_GREEN=$(tput setaf 2; tput setab 7)
_BLUE=$(tput setaf 4; tput setab 7)
_RESET=$(tput sgr0)
_BOLD=$(tput bold)

# Branch in Colab URLs:
old_url_pattern="api-examples/blob/${old_branch_name}/examples"
new_url_pattern="api-examples/blob/${new_branch_name}/examples"

styled_old_url="api-examples/blob/${_RED}${old_branch_name}${_RESET}/examples"
styled_new_url="api-examples/blob/${_GREEN}${new_branch_name}${_RESET}/examples"

# Git branches:
old_git_branch_pattern='GIT_BRANCH=\\"'"${old_branch_name}"'\\"'
new_git_branch_pattern='GIT_BRANCH=\\"'"${new_branch_name}"'\\"'

styled_old_git_branch='GIT_BRANCH=\\"'"${_RED}${old_branch_name}${_RESET}"'\\"'
styled_new_git_branch='GIT_BRANCH=\\"'"${_GREEN}${new_branch_name}${_RESET}"'\\"'

echo -e "\nFix branch in Colab URLs..."

files=$(find . -name "*.ipynb" -exec grep -H "$old_url_pattern" {} \; | cut -d: -f1)
for file in $files; do
    echo -e "File $file:"
    echo -e "    Updating URLs: ${styled_old_url} -> ${styled_new_url}"
    sed -i "s;$old_url_pattern;$new_url_pattern;g" $file
done

echo -e "\nFix git branches..."

files=$(find . -name "*.ipynb" -exec grep -H $old_git_branch_pattern {} \; | cut -d: -f1)
for file in $files; do
    echo -e "File $file:"
    echo -e "    Updating Git branch: ${styled_old_git_branch} -> ${styled_new_git_branch}"
    sed -i "s;$old_git_branch_pattern;$new_git_branch_pattern;g" $file
done


exit 0
