// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Auto-loaded by ROOT at startup (from the project root, where the macros are
// run). Silences the benign "-Wincomplete-umbrella" warnings that cling emits
// on macOS when a macro/header includes a libc++ C-compatibility wrapper such
// as <cmath> or <cctype>: ROOT 6.22's clang module map for libc++ is missing a
// few submodules (std.cmath, std.cctype, ...), which is harmless for the
// interpreter. The pragma applies to the cling session, so all subsequently
// interpreted #includes are covered.
{
   gInterpreter->ProcessLine(
       "#pragma clang diagnostic ignored \"-Wincomplete-umbrella\"");
}
