using Coverage
using Pkg

Pkg.test(; coverage=true)

coverage = process_folder("src")
LCOV.writefile("lcov.info", coverage)
clean_folder("src")
clean_folder("test")