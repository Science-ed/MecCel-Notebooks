push:
	@echo "Commiting..."
	@-git commit -m "New"
	@echo "Pushing..."
	@git push origin master

pull:
	@echo "Pulling..."
	@git reset --hard HEAD  
	@git pull origin $(BRANCH)