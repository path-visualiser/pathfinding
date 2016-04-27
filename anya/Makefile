ANYA_SRC = $(wildcard *java)
ANYA_CLASS = $(ANYA_SRC:.java=.class)

.PHONY: all
all: $(ANYA_CLASS)

.PHONY: clean
clean: 
	-rm -rf org
	-rm -rf pgraph
	-rm ScenarioRunner.class

$(ANYA_CLASS):
	javac -d ./ $(@:.class=.java) $(ANYA_SRC)

	
