
// Mini module that works similar:
// https://docs.python.org/3/library/argparse.html

pub const ARGDEF_REQUIRED: u32 = (1 << 0);
pub const ARGDEF_VARARGS: u32 =  (1 << 1);

pub const ARGDEF_DEFAULT: u32 = 0;

/// Argument Group, currently used for help message,
/// could be used to selectively parse arguments too.
///
/// Opaque data representing a group of arguments.
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct ArgGroup(usize);

/// Internal argument group data
struct ArgGroupData {
    name: &'static str,
    descr: &'static str,
}

struct ArgumentDef<T> {
    /// See: `ArgumentParser.add_argument` docs for descriptions of these members.
    id_short: &'static str,
    id_long: &'static str,
    metavar: &'static str,
    descr: &'static str,
    callback: Box<FnMut(&mut T, &[String]) -> Result<usize, String>>,
    nparams: usize,

    flag: u32,
    group: Option<ArgGroup>,
}

pub struct ArgumentParser<'a, T: 'a> {
    arg_handlers: Vec<ArgumentDef<T>>,
    arg_groups: Vec<ArgGroupData>,
    descr: &'static str,

    /// Generic data that Argument callbacks can store their data in.
    pub dest_data: &'a mut T,
}

impl <'a, T> ArgumentParser<'a, T> {


    pub fn add_argument_group(
        &mut self,
        name: &'static str,
        descr: &'static str,
    ) -> ArgGroup {
        let arg_group_handle = ArgGroup(self.arg_groups.len());
        self.arg_groups.push(
            ArgGroupData {
                name: name,
                descr: descr,
            }
        );
        return arg_group_handle;
    }

    /// Add a new argument definition.
    ///
    /// * `id_short` - Short flag name (single dash).
    /// * `id_long` - Long flag name (double dash).
    ///
    /// Note that Can be empty strings, but one must be non-empty.
    ///
    /// * `metavar` - Argument description for 'print_help',
    ///   typical usage includes NUMBER/FILE/DIR... etc.
    ///   Can also be blank.
    /// * `descr` - General description for 'print_help'.
    /// * `callback` - Takes the 'dest_data' arg from `ArgumentParser` and
    ///   a slice of arguments `nparams` long.
    ///   Returns the number of arguments used (must not exceed `nparams`)
    ///   or Err(string), when there is an error parsing.
    /// * `nparams`- Number of params, ignored when `ARGDEF_VARARGS` flag is set.
    ///   Used to check if there are enough parameters in advance.
    /// * `flag` - `ARGDEF_*` flags.
    pub fn add_argument(
        &mut self,
        id_short: &'static str,
        id_long: &'static str,
        descr: &'static str,
        metavar: &'static str,
        callback: Box<FnMut(&mut T, &[String]) -> Result<usize, String>>,
        nparams: usize,
        flag: u32,
        group: Option<ArgGroup>,
    )
    {
        if cfg!(debug_assertions) {
            if id_short.len() == 0 && id_long.len() == 0 {
                panic!("Both long and short arguments are empty");
            }
            if (flag & ARGDEF_VARARGS) != 0 && nparams != 0 {
                panic!("nparams argument set but ignored");
            }
            if let Some(group_id) = group {
                if group_id.0 >= self.arg_groups.len() {
                    panic!("Invalid group argument");
                }
            }
        }

        self.arg_handlers.push(
            ArgumentDef::<T> {
                id_short: id_short,
                id_long: id_long,
                metavar: metavar,
                descr: descr,
                callback: callback,
                nparams: nparams,
                flag: flag,
                group: group,
            }
            );
    }

    fn arg_handler_search(
        &self,
        arg: &String,
    ) -> Option<usize> {

        for (i, arg_handler) in (&self.arg_handlers).iter().enumerate() {
            if arg_handler.id_short == arg ||
               arg_handler.id_long == arg
            {
                return Some(i);
            }
        }

        return None;
    }

    pub fn parse(
        &mut self,
        args: &[String],
    ) -> Result<(), String> {
        let mut arg_handlers_used = vec![false; self.arg_handlers.len()];

        let mut i: usize = 0;
        while i < args.len() {
            if let Some(arg_handler_index) = self.arg_handler_search(&args[i]) {
                let arg_handler = &mut self.arg_handlers[arg_handler_index];
                arg_handlers_used[arg_handler_index] = true;

                i += 1;
                let args_for_handler;
                if (arg_handler.flag & ARGDEF_VARARGS) != 0 {
                    args_for_handler = &args[i..];
                } else {
                    args_for_handler = &args[i..(i + arg_handler.nparams)];
                    if args_for_handler.len() != arg_handler.nparams {
                        return Err(format!(
                            "Error '{}' expected {} parameters, received {}!",
                            args[i - 1], arg_handler.nparams, args_for_handler.len(),
                            ));
                    }
                }

                match (arg_handler.callback)(
                    &mut self.dest_data,
                    args_for_handler,
                ) {
                    Ok(nparams_used) => {
                        debug_assert!((arg_handler.flag & ARGDEF_VARARGS) != 0 ||
                                      arg_handler.nparams as usize >= nparams_used);
                        i += nparams_used;
                    }
                    Err(e) => {
                        return Err(format!(
                            "Error handling '{}': {}",
                            args[i - 1], e)
                        );
                    }
                }
            } else {
                return Err(format!(
                    "Error: '{}' unknown parameter!",
                    args[i],
                    ));
            }
        }

        // check all required args were used
        for (i, arg_handler) in (&self.arg_handlers).iter().enumerate() {
            if (arg_handler.flag & ARGDEF_REQUIRED) != 0 &&
               (arg_handlers_used[i] == false)
            {
                return Err(format!(
                    "Error: '{}{}{}' required argument not given!",
                    arg_handler.id_short,
                    if (arg_handler.id_short.len() != 0) &&
                       (arg_handler.id_long.len() != 0) {
                        "/"
                    } else {
                        ""
                    },
                    arg_handler.id_long,
                    ));
            }
        }

        Ok(())
    }

    pub fn print_help(&self) {
        println!("{}\n", self.descr);

        let arg_group_indices =
            if self.arg_groups.len() == 0 {
                vec![vec![0; self.arg_handlers.len()]]
            } else {
                let mut arg_group_indices: Vec<Vec<usize>> = vec![vec![]; self.arg_groups.len() + 1];
                for (i, arg_handler) in (&self.arg_handlers).iter().enumerate() {
                    let index = {
                        if let Some(group_id) = arg_handler.group {
                            group_id.0 + 1
                        } else {
                            0
                        }
                    };

                    arg_group_indices[index].push(i);
                }
                arg_group_indices
            };

        for (i, arg_indices) in arg_group_indices.iter().enumerate() {
            if i == 0 {
                println!("Options:");
            } else {
                let arg_group = &self.arg_groups[i - 1];
                println!("{}:\n", arg_group.name);
                if arg_group.descr.len() != 0 {
                    println!("    {}\n", arg_group.descr);
                }
            }
            self.print_help_arg_indices(&arg_indices);
            println!("\n");
        }
    }

    fn print_help_arg_indices(&self, arg_indices: &[usize]) {

        let mut options_max_len = 0;
        let mut options = vec![];
        for i in arg_indices {
            let arg_handler = &self.arg_handlers[*i];
            let option_str = format!(
                "{}{}{} {}",
                arg_handler.id_short,
                if (arg_handler.id_short.len() != 0) &&
                   (arg_handler.id_long.len() != 0) {
                    ", "
                } else {
                    ""
                },
                arg_handler.id_long,
                arg_handler.metavar,
                );

            if options_max_len < option_str.len() {
                options_max_len = option_str.len();
            }

            options.push(option_str);
        }

        // align, not optimal but doesn't really matter here,
        // could use string formatting?
        for option_str in &mut options {
            while option_str.len() < options_max_len {
                option_str.push(' ');
            }
        }

        // for i in arg_indices 
        for (i, option_str) in (arg_indices).iter().zip(options) {
            let arg_handler = &self.arg_handlers[*i];
            println!(
                "    {}  {}",
                option_str,
                arg_handler.descr,
                );
        }
    }
}

pub fn new<'a, T>(
    dest_data: &'a mut T,
    descr: &'static str,
) -> ArgumentParser<'a, T> {
    ArgumentParser::<T> {
        arg_handlers: vec![],
        arg_groups: vec![],
        descr: descr,
        dest_data: dest_data,
    }
}

