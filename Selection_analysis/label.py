#!/sw/bin/python3
"""Tool for marking newick-format-tree branches."""
import argparse
import sys
import re

__author__ = "Bogdan Kirilenko, 2018."


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def die(msg, rc=1):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("tree_file", type=str, help="Tree file for manipulations.")
    app.add_argument("branches", type=str, default="", help="Comma-separated list of branches to mark.")
    app.add_argument("--label", "-l", type=str, default="Foreground", help="Label to put to branch.")
    app.add_argument("--show_branches", "-s", action="store_true", dest="show_branches",
                     help="Show the possible branches.")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def get_branches(tree):
    """Return a list of branches."""
    # remove all special symbols
    filtered_str = re.sub(r'[^\w]', ' ', tree)
    # get rid of numbers
    filtered = [x for x in filtered_str.split() if not x.isdigit()]
    return filtered


def main():
    """Entry point."""
    args = parse_args()
    # read the file
    with open(args.tree_file, "r") as f:
        tree_content = f.read()
    # get branches possible
    all_branches = get_branches(tree_content)

    # if a user needs only the list of branches
    if args.show_branches or args.branches == "":
        # sys.stdout.write("Branches in {} are:\n".format(args.tree_file))
        sys.stdout.write("  ".join(sorted(all_branches)) + "\n")
        sys.exit(0)

    # get list of branches to mark
    req_branches = [x for x in args.branches.split(",") if x != ""]
    label = "{" + args.label + "}"
    for branch in req_branches:
        # if - or _ in tree: simple way
        if "_" in branch or "-" in branch:  # unique name
            # line turTru2_balAcu1
            tree_content = tree_content.replace(branch, branch + label)
        else:  # more complicated case, for example myoLuc2
            # it might be myoLuc2 as is as well as myoLuc2_myoDav1
            # there are two options
            # (myoLuc2: and ,myoLuc2:
            braced = "({}:".format(branch)
            commed = ",{}:".format(branch)
            tree_content = tree_content.replace(braced, "(" + branch + label + ":")
            tree_content = tree_content.replace(commed, "," + branch + label + ":")

    sys.stdout.write(tree_content)
    sys.exit(0)


if __name__ == "__main__":
    main()
