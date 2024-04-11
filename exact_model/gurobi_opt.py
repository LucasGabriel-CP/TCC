from optimization import Optmizer
import sys


def main(argv):
    model = Optmizer()

    instances = [f"inst{i}.txt" for i in range(31)]

    if len(argv) > 2:
        model = Optmizer()
        model.add_data(instances[int(argv[2])])
        model.lesgo(time_limit=60*15, problem=argv[1].upper())
        with open('./results.txt', 'a') as fp:
            fp.write(f'Fitness: {model.Z}, gap: {model.gap}, time: {model.runtime}, bound {model.model.ObjBound}\n')
    else:
        for inst in instances[1:]:
            model = Optmizer()
            model.add_data(inst)
            model.lesgo(time_limit=60*15, problem=argv[1].upper())
            with open('./results.txt', 'a') as fp:
                fp.write(f'Fitness: {model.Z}, gap: {model.gap}, time: {model.runtime}, bound {model.model.ObjBound}\n')
    

if __name__ == "__main__":
    main(sys.argv)

