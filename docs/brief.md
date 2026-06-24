WP3 Take-Home — CRISPR guide design for
a cancer indication
Time-box: ~1–2 days.
What we care about: how you think. We're deliberately not handing you
a recipe or a checklist — a large part of what we're evaluating is whether you can work out what
matters here and justify it. This is not a leaderboard contest.
Context — what we're building
FinalDose is building a programmable, CRISPR-based cancer therapy. Instead of targeting a
protein, it targets the cancer's DNA. A guide RNA directs the molecule to a specific DNA
sequence — a mutation that marks the cancer — and only when the molecule finds and
binds that exact sequence does it change shape and release its payload: a kill switch that
destroys the cell. Find the mutation in a cancer cell and it dies; a healthy cell carrying the normal
(wild-type) sequence is never triggered and is left untouched.
Two consequences make this a guide-RNA problem first and foremost:
The guide is what programs the therapy — aiming it at a different cancer is a change of
guide RNA, not a new drug.
Because recognition is what pulls the trigger, wherever the guide directs recognition is
where the payload fires — on the intended mutation (the therapeutic effect) or anywhere it
shouldn't (harm to a healthy cell). Effectiveness and safety both live in the guide.
This track is about guide RNA design: choosing the guide(s) that make the therapy work and
keep it safe.
Your task
Your indication is pancreatic ductal adenocarcinoma (PDAC).
Show us how you would approach designing this therapy for PDAC — from *"we want to treat
this cancer"* all the way to a concrete, defensible guide-design recommendation. We want to
see your reasoning: the choices you make and why, the data and tools you bring to bear, and
how you decide whether your design is actually any good.
There's no single right answer and no checklist we're matching against. Treat it as: if this were
your project, what would you do, and how would you convince us it's the right call?
Deliverables
Code — reproducible, short README.
Writeup (≤2 pages) — your approach, your recommendation, the trade-offs you weighed,
and honest limitations.
Ground rules
We care which considerations you identify and how rigorously you reason about them —
not coverage of some list. Go deep where it matters; say what you'd do with more
time/data/wet-lab.
Honest limitations beat a tuned black box.
Any tools / libraries / LLMs are fine — it'd be great if you can walk us through your decision
tree.