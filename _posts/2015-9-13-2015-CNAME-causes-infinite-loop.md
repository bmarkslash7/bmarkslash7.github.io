---
layout: post
title: CNAME leads to infinite loop
---

I'm trying to get my website to work with redirects.  So far, I have the redirect from stieha.com to my github site.  I had CNAME set as stieha.com, which seemed to have caused an infinite loop as the website bounced from stieha.com to github back to stieha.com.  I need to be more systematic about this AND account for the lag between update and DNS resolution.
