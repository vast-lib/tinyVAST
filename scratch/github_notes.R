

# To fix an error in github commits
#  * Type "cmd" in Windows search-bar to open Windows powershell
#  * Change directory to location for github e.g., type "cd C:\Users\James.Thorson\Desktop\Git\tinyVAST"
#  * Find commit immediately prior to bad commit and click "Copy full SHA" to get SHA hash
#  * In powershell, type "git reset --hard [SHA hash]"
#  * Open Windows Explorer and browse local copy to confirm that it reverted to correct state
#  * Open Github for Windows and confirm that there are "Fetch" updates, but do not fetch (fetching these would revert to current commit for branch)
#  * In powershell, type "git push origin [branch] --force" to update GitHub cloud back to that local state

