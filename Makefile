.PHONY: all build test coverage clean

all: build

build:
	cargo build --release

test: build
ifeq ($(OS),Windows_NT)
	@echo "Skipping test on Windows"
else
	cargo run --release -- tests/sids.dbc tests/sids_out.dbf
	cmp tests/sids_out.dbf tests/sids.dbf
endif

coverage:
ifeq ($(OS),Windows_NT)
	@echo "Skipping coverage on Windows"
else
	@echo "=== Instalando cargo-tarpaulin ==="
	cargo install cargo-tarpaulin || true
	@echo "=== Validando meta de 97% de cobertura com tarpaulin ==="
	cargo tarpaulin --fail-under 97 --out Html
endif

clean:
	cargo clean
	rm -f tests/sids_out.dbf tarpaulin-report.html
