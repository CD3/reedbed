SHELL := /bin/bash

all: plot

output/reedbed/Tvst.txt:
	mkdir -p output/reedbed/
	powerconf render CONFIG.yml CONFIG-reedbed.json.template /dev/stdout | time cargo run --release -- start -f json	-t json > output/reedbed/Tvst.txt 2> reedbed.stderr

output/retina-therm/Tvst-CW.txt: CONFIG.yml CONFIG-retina-therm-CW.yml.template
	mkdir -p output/retina-therm/
	powerconf render CONFIG.yml CONFIG-retina-therm-CW.yml.template CONFIG-retina-therm-CW.yml
	retina-therm temperature-rise CONFIG-retina-therm-CW.yml

output/retina-therm/Tvst-MP.txt: output/retina-therm/Tvst-CW.txt CONFIG-retina-therm-MP.yml.template
	powerconf render CONFIG.yml CONFIG-retina-therm-MP.yml.template CONFIG-retina-therm-MP.yml
	retina-therm multiple-pulse CONFIG-retina-therm-MP.yml

.PHONY: check
check:
	@which powerconf >/dev/null || echo "powerconf cli not found. Please install (e.g. `uv tool install powerconf`) and/or make sure the directory with the executable is in your PATH"
	@which retina-therm >/dev/null || echo "retina-therml cli not found. Please install (e.g. `uv tool install powerconf`) and/or make sure the directory with the executable is in your PATH"
	@which gnuplot >/dev/null || echo "gnuplot not found. Please install."

.PHONY: run
run: output/reedbed/Tvst.txt output/retina-therm/Tvst-MP.txt


.PHONY: clean
clean:
	rm -rf output *.stderr

.PHONY: plot
plot: run
	cd output && gnuplot --persist -e 'plot "reedbed/Tvst.txt", "retina-therm/Tvst-MP.txt" using 1:2 with line title "rise"'
