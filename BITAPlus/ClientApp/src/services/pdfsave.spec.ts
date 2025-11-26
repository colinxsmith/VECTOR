import { TestBed } from '@angular/core/testing';

import { Pdfsave } from './pdfsave';

describe('Pdfsave', () => {
  let service: Pdfsave;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(Pdfsave);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
