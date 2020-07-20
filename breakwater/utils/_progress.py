import sys

class ProgressBar:
    """ Simple progress bar to show progress of a process

    Parameters
    ----------
    number : int
        width of the progress bar
    task : str
        text in front of the progress bar, for instance *Computing*
    max_width : int, optional, default: 50
        maximum width of the progress bar
    """

    def __init__(self, number, task, max_width=50):
        """ See help(ProgressBar) for more info """
        # set attibutes
        self.total = number
        self.task = task
        self.progress = 0

        # set up the progress bar
        self.file = sys.stderr

        # private variable to track if an update should happen
        self._counter = 1
        self._current_number = 0

        # check if maximum width is exceeded by the input
        if number > max_width:
            # if exceeded set width to max allowed width
            self.step = max_width/number
            self.width = max_width
        else:
            # if not exceeded, do nothing
            self.step = 1
            self.width = number

        # print first line
        line = self._line(progress=0)
        print('\r' + line, file=self.file, end='')

    def _line(self, progress):
        """ Print progress bar """
        # compute percentage progress
        bar = '[' + '=' * progress + '.' * (self.width - progress) + ']'
        line = f'{self.task} {bar} {self._current_number}/{self.total}'
        return line

    def next(self):
        """ Add next character to progress bar """
        # create a new line
        self.progress += self.step
        self._current_number += 1

        # check if new update should occur
        if self.progress >= self._counter:
            # update progress bar
            line = self._line(progress=self._counter)
            print('\r' + line, file=self.file, end='')
            self._counter += 1
        else:
            # self._counter-1 because only fraction next to bar is updated
            line = self._line(progress=self._counter-1)
            print('\r' + line, file=self.file, end='')

    def finish(self):
        """ Finish the progress bar """
        # print last part of progress bar
        line = self._line(progress=self._counter)
        print('\r' + line, file=self.file, end='')

        # end progress bar by ending the line
        self.file.write("\n")
