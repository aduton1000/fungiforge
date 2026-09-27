"""The site config must repoint every process label that uses the base image.

conf/base.config gives those labels `container = params.fungiforge_image`, a Docker Hub
reference. A site install builds that image locally as a .sif and overrides the labels in
share/site.config.example. A label present in base.config but missing from that selector
therefore keeps the Hub reference and is pulled mid-run — which fails outright when the
repository is private. Stage 16 (label `cohort`) did exactly that after hours of work.
"""
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def base_image_labels():
    text = open(os.path.join(REPO, "conf", "base.config")).read()
    return [m.group(1) for m in
            re.finditer(r"^\s*withLabel:\s*([a-z_]+).*params\.fungiforge_image", text, re.M)]


def site_config_selector():
    text = open(os.path.join(REPO, "share", "site.config.example")).read()
    m = re.search(r"withLabel:\s*'([a-z_|]+)'\s*\{\s*container[^}]*fungiforge-[0-9]", text)
    assert m, "no base-image label selector found in share/site.config.example"
    return m.group(1).split("|")


def test_base_config_has_base_image_labels():
    assert len(base_image_labels()) >= 5


def test_site_config_covers_every_base_image_label():
    missing = sorted(set(base_image_labels()) - set(site_config_selector()))
    assert not missing, (
        "share/site.config.example does not repoint these base-image labels: %s. Each would be "
        "pulled from Docker Hub during a run instead of using the locally built .sif."
        % ", ".join(missing))


def test_site_config_selector_has_no_unknown_labels():
    extra = sorted(set(site_config_selector()) - set(base_image_labels()))
    assert not extra, ("share/site.config.example repoints labels that do not use the base image "
                       "in conf/base.config: %s" % ", ".join(extra))
