# quick testing -----------------------------------------------------------
pkgdown::clean_site(pkg = ".")
pkgdown::init_site(pkg = ".")
pkgdown::build_home_index()
pkgdown::preview_page("index.html")
pkgdown::build_article(name = "HowTo")
pkgdown::preview_page("articles/HowTo.html")

# cleanup start -----------------------------------------------------------
pkgdown::clean_site(pkg = ".")
pkgdown::init_site(pkg = ".")

# index -------------------------------------------------------------------
pkgdown::build_home(preview = TRUE)
pkgdown::build_news(preview = TRUE)


# reference ---------------------------------------------------------------
# source("pkgdown/02-pkgdown-add-to-yalm-reference.R")
pkgdown::build_reference_index()
pkgdown::build_reference()
pkgdown::preview_site(path = "/reference")


# rticles -----------------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
pkgdown::build_article("Test")
pkgdown::build_articles_index()
pkgdown::build_articles()
pkgdown::preview_site(path = "/articles")


# build -------------------------------------------------------------------
pkgdown::build_site(install=F)

pkgdown::deploy_to_branch(lazy = TRUE, clean = FALSE)




##### NEW PKGDOWN WORKFLOW FOR DEPLOYING WITHOUT REBUILDING ARTICLES LOCALLY ####

# in a terminal in the root of the project run
# git fetch origin
# git worktree add docs gh-pages
#This will now update the files inside 'docs/'
# It will skip articles because it sees the existing HTML files
pkgdown::build_site(lazy = TRUE)

deploy_lazy <- function(message = "Update site") {
  # Check if docs folder exists
  if (!dir.exists("docs")) stop("The 'docs' folder is missing. Did you set up the git worktree?")
  
  message("Building site (Lazy)...")
  pkgdown::build_site(lazy = TRUE)
  
  message("Deploying to gh-pages...")
  # The -C flag runs the git command inside the 'docs' directory
  system("git -C docs add .")
  
  # We suppress warnings here in case there is 'nothing to commit'
  suppressWarnings(
    system(sprintf('git -C docs commit -m "%s"', message))
  )
  
  system("git -C docs push origin gh-pages")
  message("Deployment Complete!")
}

deploy_lazy()
