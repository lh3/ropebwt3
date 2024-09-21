#!/usr/bin/env k8

const rb3_version = "3.7-r226";

/**************
 * From k8.js *
 **************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

/*************
 * Utilities *
 *************/

function rb3_cmd_mapflt(args)
{
	let opt = { max_diff:5, gap_size:50 };
	for (const o of getopt(args, "d:g:", [])) {
		if (o.opt == 'd') opt.max_diff = parseInt(o.arg);
		else if (o.opt == 'g') opt.gap_size = parseInt(o.arg);
	}
	if (args.length < 2) {
		print("Usage: rb3tools.js mapflt [options] <maxHap> <in.e2e>");
		print("Options:");
		print(`  -d INT      max edit distance [${opt.maxdiff}]`);
		print(`  -g INT      close a gap up to INT [${opt.gap_size}]`);
		return 1;
	}
	const max_hap = parseInt(args[0]);
	let ctg0 = "", st0 = 0, en0 = 0, gap = 0;
	let ctg1 = "", st1 = 0, en1 = 0, n_hap = 0;
	for (const line of k8_readline(args[1])) {
		let m;
		if ((m = /^QS\t(\S+):(\d+)-(\d+)\t/.exec(line)) != null) {
			ctg1 = m[1], st1 = parseInt(m[2]) - 1, en1 = parseInt(m[3]), n_hap = 0;
		} else if ((m = /^QH\t(\d+)\t(\d+)\t(\d+)/.exec(line)) != null) {
			if (n_hap > max_hap) continue;
			if (parseInt(m[3]) <= opt.max_diff)
				n_hap += parseInt(m[1]);
		} else if (line === "//") {
			if (n_hap > 0 && n_hap <= max_hap) continue;
			if (ctg1 != ctg0 || st1 > en0 + opt.gap_size) {
				if (ctg0 != "") print(ctg0, st0, en0, gap);
				ctg0 = ctg1, st0 = st1, en0 = en1, gap = 0;
			} else {
				gap += st1 > en0? st1 - en0 : 0;
				en0 = en0 > en1? en0 : en1;
			}
		}
	}
	if (ctg0 != "") print(ctg0, st0, en0, gap);
}

function rb3_cmd_call(args)
{
	let opt = { };
	let re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g;
	for (const o of getopt(args, "", [])) {
	}
	if (args.length < 2) {
		print("Usage: rb3tools.js call [options] <nHap> <in.e2e>");
		print("Options:");
		return 1;
	}
	const max_hap = parseInt(args[0]);

	class Allele {
		constructor(cnt, score, ed) {
			this.cnt = cnt, this.score = score, this.ed = ed, this.acc = 0; // acc is the number of haplotypes up to score
		}
	}

	class LocalVar {
		constructor(st, en, aid, ref, alt) {
			this.st = st, this.en = en, this.aid = aid, this.ref = ref, this.alt = alt;
			this.key = `${this.st}-${this.ref}-${this.alt}`;
		}
	}

	class Variant {
		constructor(ctg, off, len, w) {
			this.ctg = ctg, this.st = off + w.st, this.en = off + w.en, this.ref = w.ref, this.alt = w.alt;
			this.end_dist = w.st > len - w.en? w.st : len - w.en;
			this.key = `${this.ctg}-${this.st}-${this.ref}-${this.alt}`;
			this.ac_real = this.ac_ambi = this.ac_flt = this.flt_score_gap = 0;
		}
		toString() {
			let info = [`AC_GOOD=${this.ac_real}`, `AC_AMBI=${this.ac_ambi}`, `AC_FLT=${this.ac_flt}`, `FLT_GAP=${this.flt_score_gap}`];
			let flt = this.ac_real > 0 && this.ac_flt == 0? "PASS" : "DUP";
			let ref, alt, pos;
			if (this.ref.length == this.alt.length) { // SNP
				pos = this.st + 1, ref = this.ref, alt = this.alt;
			} else {
				pos = this.st, ref = `N${this.ref}`, alt = `N${this.alt}`;
			}
			let t = [this.ctg, pos, ".", ref, alt, 60, flt, info.join(";")];
			return t.join("\t");
		}
	}

	let v = [], a = [], al = [], ctg1 = "", st1 = 0, en1 = 0;
	for (const line of k8_readline(args[1])) {
		let m;
		if ((m = /^QS\t(\S+):(\d+)-(\d+)\t/.exec(line)) != null) {
			ctg1 = m[1], st1 = parseInt(m[2]) - 1, en1 = parseInt(m[3]), a = [], al = [];
		} else if ((m = /^QH\t(\d+)\t(\d+)\t(\d+)\t(\S+)/.exec(line)) != null) {
			const cnt = parseInt(m[1]), score = parseInt(m[2]), ed = parseInt(m[3]), cs = m[4];
			let x = 0;
			while ((m = re_cs.exec(cs)) != null) {
				if (m[1] == ':') {
					x += parseInt(m[2]);
				} else if (m[1] == '*') {
					a.push(new LocalVar(x, x + 1,   al.length, m[2][0].toUpperCase(), m[2][1].toUpperCase()));
					++x;
				} else if (m[1] == '+') {
					const len = m[2].length;
					a.push(new LocalVar(x, x + len, al.length, m[2].toUpperCase(), ""));
					x += len;
				} else if (m[1] == '-') {
					a.push(new LocalVar(x, x,       al.length, "", m[2].toUpperCase()));
				}
			}
			al.push(new Allele(cnt, score, ed));
		} else if (line === "//") {
			let n_hap = 0;
			// calculate al[].acc; this assumes al[] is sorted by al[].score
			for (let i = 1, j = 0; i <= al.length; ++i) {
				if (i == al.length || al[i].score != al[j].score) {
					for (let k = j; k < i; ++k)
						n_hap += al[k].cnt;
					for (let k = j; k < i; ++k)
						al[k].acc = n_hap;
					j = i;
				}
			}
			// score_cutoff is the score of the max_hap-th haplotype
			let score_cutoff = 0;
			for (let i = 0; i < al.length; ++i) {
				if (al[i].acc >= max_hap) {
					score_cutoff = al[i].score;
					break;
				}
			}
			// merge calls
			let b = [];
			a.sort(function(x, y) { return x.st != y.st? x.st - y.st : x.key == y.key? 0 : x.key < y.key? -1 : 1; });
			for (let i = 1, j = 0; i <= a.length; ++i) {
				if (i == a.length || a[i].key != a[j].key) {
					let u = new Variant(ctg1, st1, en1 - st1, a[j]), max_sc = 0;
					for (let k = j; k < i; ++k) {
						const t = al[a[k].aid];
						if (t.score > score_cutoff || (t.score == score_cutoff && t.acc <= max_hap)) {
							u.ac_real += t.cnt;
						} else if (t.score == score_cutoff) {
							max_sc = t.score;
							u.ac_ambi += t.cnt;
						} else {
							max_sc = max_sc > t.score? max_sc : t.score;
							u.ac_flt += t.cnt;
						}
					}
					u.flt_score_gap = score_cutoff - max_sc;
					print(u);
					b.push(u);
					j = i;
				}
			}
		}
	}
}

/****************
 * main functon *
 ****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: rb3tools.js <command> [arguments]");
		print("Commands:");
		print("  call           call small variants");
		print("  mapflt         generate mappability filter");
		print("  version        print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == "mapflt") rb3_cmd_mapflt(args);
	else if (cmd == "call") rb3_cmd_call(args);
	else if (cmd == "version") {
		print(rb3_version);
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
