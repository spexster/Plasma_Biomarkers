getwd()
getwd()
git remote -v

getwd()  # Should return ".../Plasma_Biomarkers"

list.files(recursive = TRUE)  # Should list all created folders

# Go back to the parent directory
setwd("..")
getwd()  # Should now be "C:/Users/farza/Documents/Plasma_Biomarkers"

# Remove the duplicate folder if it exists (optional)
unlink("Plasma_Biomarkers/Plasma_Biomarkers", recursive = TRUE)

# Create directories in the correct location
folders <- c(
  "data/raw", "data/processed", "data/external",
  "R", "notebooks", "output/figures", "output/tables",
  "docs", "renv"
)

for (folder in folders) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

list.files(recursive = TRUE)
# You should now see all your folders listed

getwd()

setwd("C:/Users/farza/Documents/Plasma_Biomarkers")

list.dirs()

getwd()

list.dirs(path = ".", full.names = TRUE, recursive = FALSE)

# Set working directory to your project
setwd("C:/Users/farza/Documents/Plasma_Biomarkers")

# Remove any problematic folders that might exist
unlink("Plasma_Biomarkers", recursive = TRUE)  # Remove duplicate if exists

# Define folders to create
folders <- c(
  "data/raw", 
  "data/processed", 
  "data/external",
  "R", 
  "notebooks", 
  "output/figures", 
  "output/tables",
  "docs", 
  "renv"
)

# Create directories with full paths
for (folder in folders) {
  dir.create(file.path(getwd(), folder), 
             recursive = TRUE, 
             showWarnings = FALSE)
}

# Define folders to create
folders <- c(
  "data/raw", 
  "data/processed", 
  "data/external",
  "R", 
  "notebooks", 
  "output/figures", 
  "output/tables",
  "docs", 
  "renv"
)

# Create directories with full paths
for (folder in folders) {
  dir.create(file.path(getwd(), folder), 
             recursive = TRUE, 
             showWarnings = FALSE)
}

r
# Check all directories
list.dirs(recursive = FALSE)  # Should show top-level folders
list.dirs(recursive = TRUE)   # Should show full structure

r
setwd("C:/Users/farza/Documents/Plasma_Biomarkers")
dir.create("test_folder")
list.files()  # Verify test_folder appears

unlink("test_folder", recursive = TRUE)

r
# Initialize renv (after installing it)
renv::init()

# Try committing again
system('git add .')
system('git commit -m "Initial folder structure"')

r
# Type 'y' and press Enter when prompted
renv::init()

r
# Add all files to Git
system("git add .")

# Make your first commit
system('git commit -m "Initial project setup with folder structure"')

# Verify the commit
system("git log")  # Should show your commit

# Check all folders exist
all(
  dir.exists("data/raw"),
  dir.exists("R"),
  dir.exists("output/figures"),
  file.exists(".gitignore"),
  file.exists("renv.lock")
)  # Should return TRUE

# Check Git status

# 1. First verify the files were properly staged for removal
system('git status')  # Should show .Rproj.user/ as "deleted" in changes to be committed

# 2. Commit the change with a clear message
system('git commit -m "Remove RStudio temporary files from version control"')

# 3. Verify they're no longer tracked
system('git status')  # Should no longer show the .Rproj.user/ modifications

writeLines(c(".Rhistory", ".RData", ".Rproj.user/"), ".gitignore", append = TRUE)

# Read existing .gitignore (if any)
current_ignore <- if(file.exists(".gitignore")) readLines(".gitignore") else character(0)

# Add new patterns (avoiding duplicates)
new_patterns <- c(".Rhistory", ".RData", ".Rproj.user/")
updated_ignore <- unique(c(current_ignore, new_patterns))

# Write back to file
writeLines(updated_ignore, ".gitignore")

# Check .gitignore content
file.show(".gitignore")

# Double-check Git status
system("git status")

# See if .gitignore needs to be committed
system("git diff -- .gitignore")

r
# 1. Check tracked files (should show no RStudio temp files)
system('git ls-files | grep -E "\.Rproj|\.Rhistory|\.RData"')  # Should return nothing

# 2. Check ignore rules are effective
system('git check-ignore -v .Rproj.user/some_file')  # Should match your ignore pattern

# Stage and commit .gitignore explicitly (though not strictly needed)
system('git add .gitignore')
system('git commit -m "Confirm RStudio files are ignored"')

r
# 1. List all tracked files (should show no RStudio temp files)
system('git ls-files')

# 2. Create a test file to confirm ignore rules work
write("test", ".Rproj.user/testfile.txt")
system('git status')  # Should NOT show the file
unlink(".Rproj.user/testfile.txt")  # Clean up