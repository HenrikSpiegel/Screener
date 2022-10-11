import argparse
from functions import submit2



default_args = dict(
        modules = "tools anaconda3/2021.11",
        runtime = 60,
        cores = 10,
        ram=40,
        group="dtu_00009",
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--command", required=True, help="path to script")
    parser.add_argument("--name", default="add_to_que_default", help="jobname")
    parser.add_argument("--test", action="store_true", help="Prints the qsub job-script to stdout and doesn't add_to_que")
    parser.add_argument("-o", nargs="+", help="options in kw:arg kwargs passed to submit2")
    args = parser.parse_args()


    final_args = default_args.copy()
    user_args = {
        "jobname": args.name,
        'command': args.command
    }
    
    if args.o:
        for kwarg in args.o:
            kw, arg = kwarg.split(":")
            user_args.update({kw: arg})
    final_args.update(user_args)
    out = submit2(**final_args, test=args.test)
    print(out)


    # jobname=jobname,
    # output = os.path.join("logs", jobname),
    # error = os.path.join("logs",jobname)