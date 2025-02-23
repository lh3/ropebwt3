#!/usr/bin/env k8

const rb3_version = "3.9-r259";

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
		if (o.opt == '-d') opt.max_diff = parseInt(o.arg);
		else if (o.opt == '-g') opt.gap_size = parseInt(o.arg);
	}
	if (args.length < 2) {
		print("Usage: rb3tools.js mapflt [options] <maxHap> <in.e2e>");
		print("Options:");
		print(`  -d INT      max edit distance [${opt.max_diff}]`);
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

function rb3_e2e_read1(f, buf, thres1, thres2)
{
	let r = { c1:0, c2:0, ctg:null, st:-1, en:-1 };
	while (f.readline(buf) >= 0) {
		let m, line = buf.toString();
		if ((m = /^QS\t(\S+):(\d+)-(\d+)\t/.exec(line)) != null) {
			r.ctg = m[1], r.st = parseInt(m[2]) - 1, r.en = parseInt(m[3]);
		} else if ((m = /^QH\t(\d+)\t(\d+)\t(\d+)/.exec(line)) != null) {
			let ed = parseInt(m[3]), cnt = parseInt(m[1]);
			if (ed <= thres1) r.c1 += cnt;
			if (ed <= thres2) r.c2 += cnt;
		} else if (line == "//") {
			break;
		}
	}
	return r.ctg != null? r : null;
}

function rb3_cmd_mapflt2(args)
{
	let opt = { max_rdiff:3, max_pdiff:7, gap_size:50 };
	for (const o of getopt(args, "p:r:g:", [])) {
		if (o.opt == '-p') opt.max_pdiff = parseInt(o.arg);
		else if (o.opt == '-r') opt.max_rdiff = parseInt(o.arg);
		else if (o.opt == '-g') opt.gap_size = parseInt(o.arg);
	}
	if (args.length < 3) {
		print("Usage: rb3tools.js mapflt2 [options] <maxHap> <in.ref.e2e> <in.pan.e2e>");
		print("Options:");
		print(`  -r INT      max edit distance for reference [${opt.max_rdiff}]`);
		print(`  -p INT      max edit distance for pangenome [${opt.max_pdiff}]`);
		print(`  -g INT      close a gap up to INT [${opt.gap_size}]`);
		return 1;
	}
	const max_hap = parseInt(args[0]);
	let fr = new File(args[1]), fp = new File(args[2]);
	let buf = new Bytes();
	let ctg0 = "", st0 = 0, en0 = 0, gap = 0, r;
	while ((r = rb3_e2e_read1(fr, buf, opt.max_rdiff, opt.max_pdiff)) != null) {
		let p = rb3_e2e_read1(fp, buf, opt.max_rdiff, opt.max_pdiff);
		if (p == null) throw Error("more records in the reference e2e file");
		if (r.ctg != p.ctg || r.st != p.st || r.en != p.en) throw Error("inconsistent coordinate");
		let flt = false;
		if (r.c1 == 1 && p.c1 > 0 && p.c1 <= max_hap) { // reference count is exactly 1; pangenome count no more than the threshold
			if (r.c2 == 1 && p.c2 > max_hap) flt = true;
		} else flt = true;
		if (flt) {
			if (r.ctg != ctg0 || r.st > en0 + opt.gap_size) {
				if (ctg0 != "") print(ctg0, st0, en0, gap);
				ctg0 = r.ctg, st0 = r.st, en0 = r.en, gap = 0;
			} else {
				gap += r.st > en0? r.st - en0 : 0;
				en0 = en0 > r.en? en0 : r.en;
			}
		}
	}
	if (ctg0 != "") print(ctg0, st0, en0, gap);
	buf.destroy();
	fr.close(); fp.close();
}

function rb3_cmd_call(args)
{
	let opt = { dbg:false, ambi_range:4, drop_score:12, max_gced:5, keep_supp1:false, flag_conflict:false };
	let re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g;
	for (const o of getopt(args, "r:a:d:1c", ["dbg"])) {
		if (o.opt == "--dbg") opt.dbg = true;
		else if (o.opt == "-r") opt.drop_score = parseInt(o.arg);
		else if (o.opt == "-a") opt.ambi_range = parseInt(o.arg);
		else if (o.opt == "-d") opt.max_gced = parseInt(o.arg);
		else if (o.opt == "-1") opt.keep_supp1 = true;
		else if (o.opt == "-c") opt.flag_conflict = true;
	}
	if (args.length < 2) {
		print("Usage: rb3tools.js call [options] <nHap> <in.e2e>");
		print("Options:");
		print(`  -d INT     max gap-compressed edit distance [${opt.max_gced}]`);
		print(`  -a INT     filter a variant if score within cutoff-INT [${opt.ambi_range}]`);
		print(`  -r INT     drop alignments with score lower than cutoff-INT [${opt.drop_score}]`);
		print(`  -1         keep variants supported by one k-mer only`);
		print(`  -c         output the CONFLICT filter`);
		return 1;
	}
	const max_hap = parseInt(args[0]);

	print(`##fileformat=VCFv4.2`);
	print(`##source=rb3tools-${rb3_version}`);
	print(`##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of alternate allele">`);
	print(`##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of samples">`);
	print(`##INFO=<ID=AC_AMBI,Number=A,Type=Integer,Description="Number of ambiguous alleles">`);
	print(`##INFO=<ID=AN_AMBI,Number=1,Type=Integer>`);
	print(`##INFO=<ID=AC_DUP,Number=A,Type=Integer,Description="Number of duplicate alleles">`);
	print(`##INFO=<ID=AN_DUP,Number=1,Type=Integer>`);
	print(`##INFO=<ID=RSCORE,Number=1,Type=Integer,Description="Relative k-mer alignment score">`);
	print(`##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of supporting k-mers">`);
	print(`##FILTER=<ID=LOWCONF,Description="Low confidence">`);
	print(`##FILTER=<ID=AMBI,Description="Ambiguous">`);
	print(`##FILTER=<ID=DUP,Description="Likely caused by duplications">`);
	print(`##FILTER=<ID=SUPPORT1,Description="Supported by one k-mer only">`);
	if (opt.flag_conflict)
		print(`##FILTER=<ID=CONFLICT,Description="Conflictive with a better k-mer alignment">`);
	print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO");

	class Allele {
		constructor(cnt, score, ed) {
			this.cnt = cnt, this.score = score, this.ed = ed, this.acc = 0; // acc is the number of haplotypes up to score
			this.type = -1; // -1 for unset
		}
	}

	class KmerVar {
		constructor(st, en, aid, ref, alt) {
			this.st = st, this.en = en, this.aid = aid, this.ref = ref, this.alt = alt;
			this.key = `${this.st}-${this.ref}-${this.alt}`;
		}
	}

	class Variant {
		constructor(kmer_id, ctg, off, len, w) {
			this.kmer_id = kmer_id, this.ctg = ctg, this.st = off + w.st, this.en = off + w.en, this.ref = w.ref, this.alt = w.alt;
			this.end_dist = w.st < len - w.en? w.st : len - w.en;
			this.conflict_flt = false;
			this.key = `${this.ctg}-${this.st}-${this.ref}-${this.alt}`;
			this.ac_real = this.ac_ambi = this.ac_flt = 0;
			this.an_real = this.an_ambi = this.an_flt = 0;
			this.rel_score = 0;
			this.n_support = 1;
			this.type = -1;
		}
		toString() {
			let info = [`AC=${this.ac_real}`, `AN=${this.an_real}`, `AC_AMBI=${this.ac_ambi}`, `AN_AMBI=${this.an_ambi}`,
						`AC_DUP=${this.ac_flt}`, `AN_DUP=${this.an_flt}`, `RSCORE=${this.rel_score}`, `SUPPORT=${this.n_support}`];
			let flt = [];
			if (this.type > 0) flt.push(this.type == 1? "LOWCONF" : this.type == 2? "AMBI" : "DUP");
			if (!opt.keep_supp1 && this.n_support < 2) flt.push("SUPPORT1");
			if (opt.flag_conflict && this.conflict_flt) flt.push("CONFLICT");
			if (flt.length == 0) flt.push("PASS");
			let ref, alt, pos;
			if (this.ref.length == this.alt.length) { // SNP
				pos = this.st + 1, ref = this.ref, alt = this.alt;
			} else {
				pos = this.st, ref = `N${this.ref}`, alt = `N${this.alt}`;
			}
			let t = [this.ctg, pos, ".", ref, alt, 60, flt.join(";"), info.join(";")];
			return t.join("\t");
		}
	}

	let kmer_id = 0, vcf = [], a = [], al = [], ctg1 = "", st1 = 0, en1 = 0;
	for (const line of k8_readline(args[1])) {
		let m;
		if ((m = /^QS\t(\S+):(\d+)-(\d+)\t/.exec(line)) != null) {
			ctg1 = m[1], st1 = parseInt(m[2]) - 1, en1 = parseInt(m[3]), a = [], al = [];
		} else if ((m = /^QH\t(\d+)\t(\d+)\t(\d+)\t(\S+)/.exec(line)) != null) {
			const cnt = parseInt(m[1]), score = parseInt(m[2]), ed = parseInt(m[3]), cs = m[4];
			let x = 0, gced = 0, b = [];
			while ((m = re_cs.exec(cs)) != null) {
				if (m[1] == ':') {
					x += parseInt(m[2]);
				} else if (m[1] == '*') {
					b.push(new KmerVar(x, x + 1,   al.length, m[2][0].toUpperCase(), m[2][1].toUpperCase()));
					++x, ++gced;
				} else if (m[1] == '+') {
					const len = m[2].length;
					b.push(new KmerVar(x, x + len, al.length, m[2].toUpperCase(), ""));
					x += len, ++gced;
				} else if (m[1] == '-') {
					b.push(new KmerVar(x, x,       al.length, "", m[2].toUpperCase()));
					++gced;
				}
			}
			if (gced <= opt.max_gced) {
				for (let i = 0; i < b.length; ++i)
					a.push(b[i]);
				al.push(new Allele(cnt, score, ed));
			}
		} else if (line === "//") {
			if (opt.dbg) print("X1", `${ctg1}:${st1+1}-${en1}`);
			// output variants that have moved out of the current window
			while (vcf.length && (vcf[0].ctg != ctg1 || vcf[0].en <= st1))
				print(vcf.shift());
			// calculate al[].acc; this assumes al[] is sorted by al[].score
			let n_hap = 0;
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
			let score_cutoff = 0, score_next = 0;
			for (let i = 0; i < al.length; ++i) {
				if (al[i].acc >= max_hap && score_cutoff == 0)
					score_cutoff = al[i].score;
				if (al[i].acc > max_hap && score_next == 0)
					score_next = al[i].score;
			}
			if (score_cutoff == 0 && al.length > 0)
				score_cutoff = al[al.length - 1].score;
			if (opt.dbg) print("X2", score_cutoff, score_next);
			// classify each allele
			let an_real = 0, an_ambi = 0, an_flt = 0;
			for (let i = 0; i < al.length; ++i) {
				let t = al[i];
				if (t.score >= score_cutoff && t.score >= score_next + opt.ambi_range) t.type = 0, an_real += t.cnt; // probably real
				else if (t.score >= score_cutoff && t.score > score_next) t.type = 1, an_real += t.cnt; // perhaps real but not reliable
				else if (t.score < score_cutoff - opt.drop_score) t.type = 4; // discard
				else if (t.score == score_next) t.type = 2, an_ambi += t.cnt; // can't tell
				else if (t.score < score_next) t.type = 3, an_flt += t.cnt; // probably false variants
			}
			an_flt += an_real + an_ambi;
			an_ambi += an_real;
			if (score_cutoff == score_next) an_real = max_hap;
			// merge calls
			a.sort(function(x, y) { return x.key == y.key? 0 : x.key < y.key? -1 : 1; });
			for (let i = 1, j = 0; i <= a.length; ++i) {
				if (i == a.length || a[j].key != a[i].key) {
					let v = new Variant(kmer_id, ctg1, st1, en1 - st1, a[j]), max_sc = 0, best_type = 4;
					for (let k = j; k < i; ++k) {
						const t = al[a[k].aid];
						best_type = best_type < t.type? best_type : t.type;
						if (t.type == 4) continue;
						else if (t.type <= 1) v.ac_real += t.cnt, v.an_real = 0;
						else if (t.type == 2) v.ac_ambi += t.cnt;
						else if (t.type == 3) v.ac_flt += t.cnt;
						max_sc = max_sc > t.score? max_sc : t.score;
					}
					if (best_type < 4) {
						v.type = best_type;
						v.rel_score = max_sc - score_cutoff;
						v.an_real = an_real, v.an_ambi = an_ambi, v.an_flt = an_flt;
						vcf.push(v);
					}
					j = i;
				}
			}
			// resolve conflicts with other k-mers
			let wcf = [];
			vcf.sort(function(x, y) { return x.st != y.st? x.st - y.st : x.key == y.key? 0 : x.key < y.key? -1 : 1; });
			for (let i = 1, j = 0; i <= vcf.length; ++i) {
				if (i == vcf.length || vcf[j].key != vcf[i].key) {
					let n_curr = 0, max_end_dist = -1, max_k = -1, n_support = 0;
					for (let k = j; k < i; ++k) {
						const v = vcf[k];
						if (v.kmer_id == kmer_id)
							++n_curr;
						if (v.end_dist > max_end_dist)
							max_end_dist = v.end_dist, max_k = k;
						n_support += v.n_support;
					}
					if (n_curr > 1 || max_k < 0) throw("Bug!");
					let v = vcf[max_k];
					v.n_support = n_support;
					if (n_curr == 0) {
						const curr_end_dist = v.st - st1 < en1 - v.en? v.st - st1 : en1 - v.en;
						if (v.end_dist < curr_end_dist)
							v.conflict_flt = true;
					}
					wcf.push(v);
					j = i;
				}
			}
			vcf = wcf;
			++kmer_id;
		}
	}
	while (vcf.length) // output the rest of calls
		print(vcf.shift());
}

function rb3_cmd_getsnp(args)
{
	let opt = { auto_only:false };
	for (const o of getopt(args, "a", [])) {
		if (o.opt == '-a') opt.auto_only = true;
	}
	if (args.length < 1) {
		print("Usage: rb3tools.js getsnp [options] <in.vcf>");
		print("Options:");
		print(`  -a         chromsome names must match /^(chr[0-9]+|[0-9]+)$/`);
		return 1;
	}
	for (const line of k8_readline(args[0])) {
		if (typeof line != "string" || line.length == 0 || line[0] == '#') continue; // header line
		let t = line.split("\t", 8);
		if (opt.auto_only && !/^(chr\d+|\d+)$/.test(t[0])) continue;
		const ref = t[3];
		let s = t[4].split(",");
		for (let i = 0; i < s.length; ++i) {
			const alt = s[i];
			if (ref.length != alt.length) continue;
			for (let k = 0; k < ref.length; ++k)
				if (ref[k] != alt[k])
					print([t[0], t[1], ref[k], alt[k]].join("-"));
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
		//print("  mapflt         generate mappability filter");
		print("  mapflt2        generate mappability filter");
		print("  getsnp         extract SNPs");
		print("  version        print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == "mapflt") rb3_cmd_mapflt(args);
	else if (cmd == "mapflt2") rb3_cmd_mapflt2(args);
	else if (cmd == "call") rb3_cmd_call(args);
	else if (cmd == "getsnp") rb3_cmd_getsnp(args);
	else if (cmd == "version") print(rb3_version);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
