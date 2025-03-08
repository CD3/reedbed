SHELL := /bin/bash

all: plot

output/: output/reedbed/finished.txt output/retina-therm/finished.txt

output/reedbed/finished.txt: CONFIG.yml CONFIG-reedbed.json.template
	powerconf run reedbed CONFIG.yml
	touch output/reedbed/finished.txt

output/retina-therm/finished.txt: CONFIG.yml CONFIG-retina-therm.yml.template
	powerconf run retina-therm CONFIG.yml
	touch output/retina-therm/finished.txt

output/pac1d/finished.txt: CONFIG.yml CONFIG-pac1d.yml.template
	powerconf run pac1d-cli CONFIG.yml
	touch output/pac1d/finished.txt

.PHONY: check-config
check-config:
	powerconf print-instances CONFIG.yml

.PHONY: check-tools
check-tools:
	@which powerconf >/dev/null || echo "powerconf cli not found. Please install (e.g. pipx install powerconf) and/or make sure the directory with the executable is in your PATH"
	@ [[ $$(powerconf version | gawk '{print $$2}' | sed '$$a0.3.0' | sort -V | head -n 1) == 0.3.0 ]] || echo "powerconf $$(powerconf version) is too old. Please install powerconf>=0.3.0"
	@which retina-therm >/dev/null || echo "retina-therm cli not found. Please install (e.g. pipx install retina-therm) and/or make sure the directory with the executable is in your PATH"
	@which awk >/dev/null || echo "awk not found. Please install."
	@which sed >/dev/null || echo "sed not found. Please install."

.PHONY: run
run: output/

.PHONY: clean
clean:
	rm -rf output config_files \*.stderr

.PHONY: plot
plot: run
	cd output && gnuplot --persist -e 'plot "reedbed/Tvst.txt", "retina-therm/Tvst.txt"
