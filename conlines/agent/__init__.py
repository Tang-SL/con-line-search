"""Line-identification agent toolkit.

A command-line toolkit (``python -m conlines.agent <step> ...``, or ``scripts/agent``)
that an AI agent (Claude Code, following .claude/skills/line-id-agent/SKILL.md) or a
person drives step by step. Every step reads and writes a *run folder*:

    run.json        inputs and state (cubes, target, redshift, aperture, settings)
    *.json          step results (inventory, context, components, features, ...)
    plots/          PNG figures for checking every step
    products/       FITS maps and tables
    notes/          the agent's own reasoning, one Markdown file per step
    log.jsonl       what ran, when, with which settings

``report`` turns the run folder into a static HTML weblog.
"""
