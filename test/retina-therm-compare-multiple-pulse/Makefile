SHELL := /bin/zsh

.PHONY: all
all: plot

.PHONY: configs
configs: CONFIG-reedbed.json CONFIG-retina-therm-cw.yml CONFIG-retina-therm-mp.yml

CONFIG-reedbed.json: CONFIG.yml
	powerconf generate CONFIG.yml -n reedbed -f json CONFIG-reedbed.json

CONFIG-retina-therm-cw.yml: CONFIG.yml
	powerconf generate CONFIG.yml -n retina-therm-cw -f yml CONFIG-retina-therm-cw.yml

CONFIG-retina-therm-mp.yml: CONFIG.yml
	powerconf generate CONFIG.yml -n retina-therm-mp -f yml CONFIG-retina-therm-mp.yml

output/reedbed/Tvst.txt: CONFIG-reedbed.json
	mkdir -p output/reedbed/
	cat CONFIG-reedbed.json | time cargo run --release -- start -f json -t json > output/reedbed/Tvst.txt 2> reedbed.stderr

output/retina-therm/Tvst-CW.txt: CONFIG-retina-therm-cw.yml
	mkdir -p output/retina-therm/
	retina-therm temperature-rise CONFIG-retina-therm-cw.yml

output/retina-therm/Tvst-MP.txt: output/retina-therm/Tvst-CW.txt CONFIG-retina-therm-mp.yml
	retina-therm multiple-pulse CONFIG-retina-therm-mp.yml

.PHONY: check
check:
	@which powerconf >/dev/null || echo "powerconf cli not found. Please install (e.g. `uv tool install powerconf`) and/or make sure the directory with the executable is in your PATH"
	@which retina-therm >/dev/null || echo "retina-therml cli not found. Please install (e.g. `uv tool install powerconf`) and/or make sure the directory with the executable is in your PATH"
	@which gnuplot >/dev/null || echo "gnuplot not found. Please install."

.PHONY: run
run: output/reedbed/Tvst.txt output/retina-therm/Tvst-MP.txt

.PHONY: clean
clean:
	rm -rf output *.stderr CONFIG-*

.PHONY: plot
plot: run
	cd output && gnuplot --persist -e 'plot "reedbed/Tvst.txt" title "reedbed", "retina-therm/MP/Tvst.txt" using 1:2 with line title "retina-therm"'
